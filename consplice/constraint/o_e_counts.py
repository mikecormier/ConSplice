from __future__ import print_function
import sys
import io
import os
import argparse 
import gzip
import copy
from cyvcf2 import VCF, Writer 
import pysam
import pandas as pd
import numpy as np
import pyranges as pr
from collections import defaultdict
import multiprocessing
import datetime
from pyfaidx import Fasta
import time
from .utils import *

#---------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
#---------------------------------------------------------------------------------------------------------------------------------

def add_o_e_counts(sub_p):


    p = sub_p.add_parser("oe-counts",
                         help = "Calculate Observed and Expected splicing variant counts for genic or intragenic regions of genes",
                         formatter_class=argparse.RawDescriptionHelpFormatter,
                         description = ("\n\t********************************************\n"
                                        "\t* ConSplice - Observed and Expected Counts *\n"
                                        "\t********************************************\n\n"
                                        "\tCalculate the Observed (O) splicing variation from genetic"
                                        "\n\tvariation in gnomAD, and calculate the Expected (E)"
                                        "\n\tsplicing variation using a matrix of Substitution"
                                        "\n\trates based on evidence of alternative splicing."  
                                        "\n\tO and E counts are based on a gene or window based region"
                                        "\n\tdesignated by the user."
                                       ),
    )

    req = p.add_argument_group("Required Arguments")
    region_req = p.add_argument_group("Required Arguments for 'Region'")

    req.add_argument(
        "--gnomad-vcf",
        metavar = "gnomAD VCF/BCF File", 
        required=True, 
        help="(Required) The path to the gnomAD vcf/bcf file. (bgzipped and tabixed required)"
    )

    req.add_argument(
        "--coverage",
        metavar = "gnomAD Coverage file",
        required = True,
        help = "(Required) The gnomAD coverage file for the gnomAD vcf file. The coverage files should be in bed format, bgzipped, and tabixed."
    )

    req.add_argument(
        "--gtf-file",
        metavar = "GTF file",
        required=True,
        help="(Required) A gtf file to get gene specific genomic coordinance from."
    )   

    req.add_argument(
        "--spliceai-vcf",
        metavar = "SpliceAI SNV prediction File", 
        required=True, 
        help="(Required) The path to SpliceAI delta score predictions for every snv (vcf file)."
    )

    req.add_argument(
        "--alt-gene-symbol",
        metavar = "Alternative Gene Symbol File",
        required = True,
        help = "(Required) A file that contains mappings between a canonical gene symbol to alternative gene symbols. NOTE: The file needs to have a header!. This script is set up to use the HGNC protein-coding gene mapping file: ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt, however, you are welcome to try others. There is no guarantee it will work with other mapping files."  
    )

    req.add_argument(
        "--fasta",
        metavar = "Fasta file",
        required=True,
        help="(Required) A fasta file to get reference alleles by position from."
    )   

    req.add_argument(
        "--seg-dups",
        metavar="Segmental Duplications File",
        required=True,
        help="(Required) The file path to the segmental duplications file. This file will be used to identify repeat regions that should be excluded from the O and E calcualtion. (We suggest using the ggd segmental duplications data package 'grch38-segmental-dups-ucsc-v1')."
    )

    req.add_argument(
        "--self-chains",
        metavar="Self Chains",
        required=True,
        help="(Required) The file path to the high identity self chain repeat file. This file will be used to identify repeat regions that should be exlucded from the O and E calcualtion. (We suggest using the ggd self chain data package 'grch38-self-chain-ucsc-v1')."
    )

    req.add_argument(
        "--out-file",
        metavar = "Output file",
        required = True,
        help = "(Required) The name of the output file to create."
    )

    req.add_argument(
        "-sm",
        "--substitution-matrix",
        metavar="SpliceAI aware Substitution Matrix",
        required=True,
        help="(Required) Path and name of the Substitution matrix  created using SpliceAI predictions and gnomAD observed varaints. (Created from sub-commmand: sub-matrix)."
    )

    req.add_argument(
        "--region-type",
        required = True,
        choices = ["gene","region"],
        help="(Required) Calcualte observed and expected counts using the entire gene as a region or intragenic regions based on a window size. (Choices = 'gene' or 'region')"
    )

    region_req.add_argument(
        "--window-size",
        metavar="Window Size",
        type = int,
        help="(Required if --region-type set to 'region') The window size in bp used for the sliding window. Region from a gene will be created using this window size"
    )

    region_req.add_argument(
        "--step-size",
        metavar="Step Size",
        type = int,
        help = "(Required if --region-type set to 'region') The step size in bp used to slide the window along a gene. The step size will be the number of bases the window is slid before a new region is created. (NOTE: To exclude overlapping windows set the step size to the same size as the window size. For example, a window size of 25 and a step size of 25 will create non-overlapping 25bp regions without missing extra bases. If the step size is larger than the window size then missing regions will arise)"
    )

    p.add_argument(
        "--chrom",
        metavar = "Chromosomes",
        nargs = "*",
        choices = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","All"],
        help = "Which chromosomes to use to build the count matrix. For autosome only choose chromosomes 1-22, for X choose only X, All for all, etc., meaning all reference chromosomes. Multiple chromosomes should be provided as each chromsome seperated by a space (chrom 1-5 exampe: 1 2 3 4 5). (Default = ALL)."
    )

    p.add_argument(
        "-cc",
        metavar = "Coverage Cutoff",
        default = 0.5,
        type=float,
        help = "A number between 0.0 and 1.0 that represents the required number of samples with coverage at each site looked at in gnomAD. Default = 0.5"
    )

    p.add_argument(
        "--coverage-label",
        metavar = "Coverage Label",
        default = "over_10",
        help = "The label/column in the coverage file to use for the selected position coverage. Default = 'over_10'"
    )

    p.add_argument(
        "--n-cpu",
        metavar = "Number of CPUs",
        default = 3,
        type = int,
        help = "The number of CPUs to use for multi-threading. (Default = 3)"
    )

    p.add_argument(
        "--repeat-score-cutoff",
        metavar = "Repeat Score Cutoff",
        default = 0.95,
        type = float,
        help = "The cutoff score to use for seg dups or self chain repeats. Regions of a gene with a seg dup or self chain at or above this cutoff will be skipped while regions of a gene below this cutoff will not be skipped during O and E score calculations. (Default = 0.95, meaning 95 percent. The score will be 0.95 for seg dups and 95.0 for self chains)"
    )

    p.add_argument(
        "--spliceai-score-type",
        choices = ["max","sum","splicing_unaware"],
        default = "sum",
        help = "(Optional) How to use the SpliceAI score. Choices = 'max', 'sum', 'splicing_unaware'. 'max' will use the max SpliceAI score for a specific variant. 'sum' will use the sum SpliceAI score for a specific variant. 'splicing_unaware' will use a single bin for SpliceAI, which is the same as an unaware splicing model. Defulat = 'sum'"
    )

    p.set_defaults(func=o_e_counts)


#---------------------------------------------------------------------------------------------------------------------------------
## Functions/Methods
#---------------------------------------------------------------------------------------------------------------------------------


def load_config_file(config_path, region_type, chrom_list, spliceai_score_type, window_size, step_size):
    """
    load_config_file
    ================
    Method to load the config values into global variants
    
    Parameters:
    -----------
    1) config_path:         (str) Path to the config file
    2) region_type:         (str) The region_type from the input arguments
    3) chrom_list:         (list) A list of chromosomes from the input arguments
    4) spliceai_score_type: (str) The score type from the input arguments
    5) window_size:         (int) The size of the window 
    6) step_size:           (int) The step size
    """

    ## Global vars
    global SCORE_TYPE
    global delta_score_bins 
    global GENOME_BUILD
    global error_log_file
    global out_log_file
    global SAI_SCORE_TYPE
    global PAR_PR


    ## Load score type
    SCORE_TYPE = spliceai_score_type
    print("\n  Using the '{}' of SpliceAI scores".format(spliceai_score_type))

    ## Load global config
    config_dict = load_config(config_path)

    ## set global variables from config
    GENOME_BUILD = config_dict["GENOME_BUILD"]
    print("\n  Genome Build set to '{}'".format(GENOME_BUILD))

    ## Create output and error log files
    file_date = datetime.datetime.now().strftime("%m-%d-%Y")
    error_log_file = (file_date +
                      ".{}".format(region_type) + 
                      ".spliceai_{}".format(spliceai_score_type) +
                      ".O_E_Counts.chr{}.".format("chr".join(chrom_list)) +
                      ".{}window".format(window_size) +
                      ".{}step".format(step_size) + 
                       config_dict["LOG_FILES"]["error_log"]
                      )

    out_log_file = (file_date +
                    ".{}".format(region_type) + 
                    ".spliceai_{}".format(spliceai_score_type) +
                    ".O_E_Counts.chr{}.".format("chr".join(chrom_list)) +
                    ".{}window".format(window_size) +
                    ".{}step".format(step_size) + 
                    config_dict["LOG_FILES"]["out_log"]
                    )
    
    print("\n  OUT LOG:   '{}'".format(out_log_file))
    print("\n  ERROR LOG: '{}'".format(error_log_file))

    ## Score bins
    if spliceai_score_type == "max":
        delta_score_bins = config_dict["SCORE_BINS"]["max_spliceai_score_bins"]
        SAI_SCORE_TYPE = "max"

    elif spliceai_score_type == "sum":
        delta_score_bins = config_dict["SCORE_BINS"]["sum_spliceai_score_bins"]  
        SAI_SCORE_TYPE = "sum"

    elif spliceai_score_type == "splicing_unaware":
        delta_score_bins = config_dict["SCORE_BINS"]["one_sum_spliceai_score_bin"]  
        SCORE_TYPE = "sum"
        SAI_SCORE_TYPE = "sum"

    print("\n  SpliceAI score bins:\n  --------------------\n\t- {}".format("\n\t- ".join(delta_score_bins)))


    PAR_PR = create_pr_par(chrom = "X", par_dict = config_dict["PAR"][GENOME_BUILD])
    
    if "X" in chrom_list:
        print("\n  PAR Regions:")
        print(PAR_PR)


def sliding_window_regions(start, end, window_size, step_size):
    """
    sliding_window_regions
    ======================
    This method will split a gene into different regions based on a sliding window and step size. 
     Each region is based on the window size. The window is slid down the gene using the step size. 
     Each step size results in a new window. For example, if the gene is ~67000 bp long, the window size
     is 1.5kb and the step size is 375 bp, then you would get ~180 overlaping regions, with each region 
     having a size of 1.5kb. The first region will start at the start site of the gene and the last
     region will end at the end site of the gene.

    Parameters:
    -----------
    1) start:           (int) The genomic start position of a gene
    2) end:             (int) The genomic end position of a gene
    3) window_size:     (int) The size of the window/region in bp to use for the sliding window
    4) step_size:       (int) The sliding step size in bp to slide the window across a gene. 

    Returns:
    ++++++++
    (list) A 2d list of regions created by sliding a window across the gene positions. Each inner list 
            have the region start pos at index 0 and the end pos as index 1
    """
    start = int(start)
    end = int(end)
    window_size = int(window_size)
    step_size = int(step_size)

    ## Start and end of first region
    ## First region will start at the start of the gene
    window_start = start
    window_end = start + (window_size - 1) ## The size of the region will include the start position to the end position. This accounts for a off by 1 error. 

    gene_regions = []

    ## Iterate over the gene range and get all regions
    while window_end < end:
        
        ## Add region
        gene_regions.append([window_start, window_end])

        ## Slide the window by the step size
        window_start += step_size
        window_end += step_size

    ## Start and end of last region
    ## Last region will end at the end of the gene 
    window_start = end - (window_size - 1) if end - (window_size - 1) > start else start
    window_end = end
    ## Add region
    gene_regions.append([window_start, window_end])

    return(gene_regions)


def get_transcript_regions(gtf_file, biotype="protein_coding"):
    """
    get_transcript_regions
    ======================
    Method to get transcript regions from a gtf file filtered by a gene/transcript biotype. 

    Parameters:
    -----------
    1) gtf_file: (str) The path to the gtf file to parse
    2) biotype:  (str) The biotype to filter transcripts for. (Default = "protein_coding") 
                        To skip biotype filtering set the biotype parameter to None

    Returns:
    ++++++++
    1) (dict) A dictionary of transcript regions. 1st key is the transcript id, 2nd layer keys include ["chrom","gene_start","gene_end","strand","gene_id","gene_name","max_exon_number", and "cds_exon_count"]  
    2) (dict) A dictionary of gene names (symbols) to gene ids 
    """
    ## Parse GTF file and create a pandas dataframe
    try:
        fh = gzip.open(gtf_file, "rt", encoding = "utf-8") if gtf_file.endswith(".gz") else io.open(gtf_file, "rt", encoding = "utf-8") 
        fh.close()
    except IOError as e:
        print("!!Error!! Unable to read the gtf file: '{}'".format(gtf_file))
        print(str(e))
        sys.exit(1)

    ## Transcript dict
    transcript_region_dict = defaultdict(lambda: defaultdict(str))
    gene_name_dict = defaultdict(lambda: defaultdict(str))

    print("\n\tGetting {} transcript info from the gtf file".format(biotype))

    ## Iterate over gtf file
    try:
        fh = gzip.open(gtf_file, "rt", encoding = "utf-8") if gtf_file.endswith(".gz") else io.open(gtf_file, "rt", encoding = "utf-8") 
    except IOError as e:
        print("!!Error!! Unable to read the gtf file: '{}'".format(gtf_file))
        print(str(e))
        sys.exit(1)

    header = ["chrom","source","feature","start","end","score","strand","frame","attribute"]
    for line in fh:

        if line[0] == "#":
            continue

        line_dict = dict(zip(header,line.strip().split("\t")))
        line_dict.update({x.strip().replace("\"","").split(" ")[0]:x.strip().replace("\"","").split(" ")[1] for x in line_dict["attribute"].strip().split(";")[:-1]})

        ## get the biotype key 
        biotype_key = "gene_biotype" if "gene_biotype" in line_dict else "gene_type" if "gene_type" in line_dict else None 

        ## Check for matching biotype
        ### If the biotype parameter is None, then don't filter by biotype
        if (biotype_key and line_dict[biotype_key] == biotype) or biotype is None:

            ## Get transcript features
            if line_dict["feature"] == "transcript":

                    transcript_region_dict[line_dict["transcript_id"]]["chrom"] = line_dict["chrom"]
                    transcript_region_dict[line_dict["transcript_id"]]["gene_start"] = int(line_dict["start"])
                    transcript_region_dict[line_dict["transcript_id"]]["gene_end"] = int(line_dict["end"])
                    transcript_region_dict[line_dict["transcript_id"]]["strand"] = line_dict["strand"]
                    transcript_region_dict[line_dict["transcript_id"]]["gene_id"] = line_dict["gene_id"]
                    transcript_region_dict[line_dict["transcript_id"]]["gene_name"] = line_dict["gene_name"]
                    transcript_region_dict[line_dict["transcript_id"]]["transcript_id"] = line_dict["transcript_id"]

            ## Get exon feature
            elif line_dict["feature"] == "exon":

                ## check for exon count
                if not transcript_region_dict[line_dict["transcript_id"]]["max_exon_number"]:
                    transcript_region_dict[line_dict["transcript_id"]]["max_exon_number"] = int(line_dict["exon_number"])

                elif int(transcript_region_dict[line_dict["transcript_id"]]["max_exon_number"]) < int(line_dict["exon_number"]):
                    transcript_region_dict[line_dict["transcript_id"]]["max_exon_number"] = int(line_dict["exon_number"])

            ## Get CDS feature
            elif line_dict["feature"] == "CDS":

                ## check for cds exon count
                if not transcript_region_dict[line_dict["transcript_id"]]["cds_exon_count"]:
                    transcript_region_dict[line_dict["transcript_id"]]["cds_exon_count"] = 0

                ## Add cds exon count 
                transcript_region_dict[line_dict["transcript_id"]]["cds_exon_count"] += 1

                    
            ## If a gene feature
            elif line_dict["feature"] == "gene":
                
                gene_name_dict[line_dict["gene_name"]]["gene_name"] = line_dict["gene_name"]
                gene_name_dict[line_dict["gene_name"]]["gene_id"] = line_dict["gene_id"]
                gene_name_dict[line_dict["gene_name"]]["strand"] = line_dict["strand"]

    fh.close()
                    
    return(transcript_region_dict, gene_name_dict)


def get_normalized_matrix_per_label(groupby_object, mutation_rate_dict):
    """
    get_normalized_matrix_per_label
    ===============================
    This method is used to get the marignal distribution for each reference allele at each bin (exp: 0.8-1.0), 
     and normalize the distriubtion using Ref to alt count divided by total from the marignal distriubtion.

    Example Matrix:
      
      D+4
              ALT
              ---
        | A | C | G | T |
       -|---|---|---|---|
       A| # | # | # | # |
    R| -|---|---|---|---|
    E| C| # | # | # | # |
    F| -|---|---|---|---|
       G| # | # | # | # |
       -|---|---|---|---|
       T| # | # | # | # | 
       ------------------

    Marignal Distribution is based on Ref allele (Row)

    Normalized marginal distribution counts will add up to 1 

    Parameters:
    -----------
    1) groupby_object:     (Pandas df) A subseted pandas dataframe that represents ref to alt counts for a specific SpliceAI delta score bin (exp: "0.0-0.2")
    2) mutation_rate_dict: (dict)      A dictionary that contains the ref allele specific mutation rates. (Updated by this function)

    Returns:
    ++++++++
    Nothing is returned. Instead, the mutation_rate_dict is dynamically updated. 
    """
    
    label = groupby_object.delta_score_bin.unique()[0]

    ## Create an emtpy matrix
    zeroton_matrix = np.zeros((len(groupby_object.ref.unique()),len(groupby_object.alt.unique())), dtype= np.float)
    zeroton_plus_singleton_matrix = np.zeros((len(groupby_object.ref.unique()),len(groupby_object.alt.unique())), dtype= np.float)
    row_index = dict(zip(groupby_object.ref.unique(),[0,1,2,3]))
    column_index = dict(zip(groupby_object.alt.unique(),[0,1,2,3]))

    ## create a matrix for ref to alt counts
    for row in groupby_object.itertuples():
        ## If ref is the same as alt (no mutation)
        if row.ref == row.alt:
            
            ## Zeroton Counts
            zeroton_matrix[row_index[row.ref], column_index[row.alt]] += row.zerotons

            ## Zeroton plus Singleton counts
            zeroton_plus_singleton_matrix[row_index[row.ref], column_index[row.alt]] += row.zerotons

        else:
            ## Zeroton Counts
            ### If ref to alt, counts = non zertons
            zeroton_matrix[row_index[row.ref], column_index[row.alt]] = row.non_zerotons

            ## Zeroton plus Singleton counts
            ### Add counts to alt for all non zerotons and non singletons
            zeroton_plus_singleton_matrix[row_index[row.ref], column_index[row.alt]] = row.non_zerotons_plus_singletons
            ### Add singleton counts to the ref to ref cell
            zeroton_plus_singleton_matrix[row_index[row.ref], column_index[row.ref]] += row.singletons

    ## Get the normalized row counts
    zeroton_row_normalized_matrix = zeroton_matrix / zeroton_matrix.sum(axis = 1)[:, np.newaxis]
    zeroton_plus_singleton_row_normalized_matrix = zeroton_plus_singleton_matrix / zeroton_plus_singleton_matrix.sum(axis = 1)[:, np.newaxis]

    ## get mutation rates based on 1 - ref to ref normalized count
    zeroton_mutation_rates = [1.0 - zeroton_row_normalized_matrix[row_index[x[0]],column_index[x[1]]] for x in zip(row_index,column_index)]
    zeroton_plus_singletons_mutation_rates = [1.0 - zeroton_plus_singleton_row_normalized_matrix[row_index[x[0]],column_index[x[1]]] for x in zip(row_index,column_index)]

    ## Mutation rate to dict
    mutation_rate_dict[label]["zeroton"] = dict(zip(groupby_object.ref.unique(),zeroton_mutation_rates))
    mutation_rate_dict[label]["zeroton_plus_singleton"] = dict(zip(groupby_object.ref.unique(),zeroton_plus_singletons_mutation_rates))


def get_by_region_o_e_counts(vcf,
                                contains_chr,
                                region_chrom,
                                region_start,
                                region_end,
                                region_strand,
                                by_pos_cov_dict,
                                coverage_cutoff,
                                vep_index,
                                spliceai_region_dict,
                                alt_symbol_dict,
                                t_id_dict,
                                transcript_id,
                                mut_rate_dict,
                                non_matching_ref_alleles,
                                non_matching_symbols,
                                ref_allele_n,
                                empty_position_info):
    """
    get_by_region_o_e_counts
    ===========================
    This method is used to get observed and expected counts for a given region. it will check for the gnomAD variants 
     between the genomic region (start and end position). Variants at any of the positions in the region that pass 
     filtering are used as an observed count for a gene at a specific bin. Expected counts are then added for each 
     position in the region based on the substitution rate for that position. 
    
    If a variant passes all filters it will be added to the splice region dictionary under the appropriate ref to alt keys, 
     and the splice region label sub-key. 
     Example of Ref to Alt key: A->T
     Exmaple splice region label: A+10

    Info added include:
        zerotons count: (AC < 1) The number of sites without a variant 
        zerotons_plus_singletons count: (AC < 2) The number of sites with either a rare singleton variant or no variant 
        singeltons count: (AC == 1) The number of sites with a singelton variant 
        non_zerotons counts: (AC > 0) The number of sites with a variant 
        non_zerotons_plus_singletons count: (AC >1) The number of sites with a variant excluding singleton variants

    Parameters:
    -----------
    1)  vcf:                (cyvcf2 obj)  A cyvcf2 object for the current vcf/bcf file being queried 
    2)  contains_chr:             (bool)  whether or not the the sequence names in the vcf start with chr or not
    3)  region_chrom:              (str)  The chromsome to look in for the position interval 
    4)  region_start:              (int)  The start position of the interval
    5)  region_end:                (int)  The end position of the interval
    6)  region_strand:             (str)  A string representing the strand the gene is on. (+ or -)
    7)  by_pos_cov_dict:          (dict)  A dictionary of var coverage for the current region
    8)  coverage_cutoff:         (float)  The coverage value at each position to use as a cutoff from the coverage file. (Between 0.0 and 1.0) 
    9)  vep_index:                (list)  A list that represents the vep csq index in the gnomad vcf 
    10) spliceai_region_dict:     (dict)  A dictionary of SpliceAI information for the current region 
    11) alt_symbol_dict:          (dict)  A dictionary of different alternative gene symbols 
    12) t_id_dict:                (dict)  A by transcript id dict that tracks the observed and expected scores for each transcript id
    13) transcript_id:             (str)  The transcript id for the current transcript being investigated
    14) mut_rate_dict             (dict)  A dictionary that contains the mutation rate information
    15) non_matching_ref_alleles:  (int)  The current count of non matching ref alleles
    16) non_matching_symobls:      (int)  The current count of non matching symbols
    17) ref_allele_n:              (int)  The current count of reference alleles that are 'N'
    18) empty_position_info:       (int)  The current count of positions with null information 


    Returns:
    ++++++++
    1) (int) Updated non_matching_ref_alleles 
    2) (int) Updated non_matching_symbols 
    3) (int) Updated ref_allele_n
    4) (int) Update empty_position_info

    
    The transcript id dict, the dictionary that tracks the observed and expected score per transcript id, 
     will be dynamicially updated without returning the dict 
    """

    def get_best_delta_bin(first_delta_bin, second_delta_bin):
        """
        get_best_delta_bin
        ==================
        Method to get the best/highest delta score bin between two bins. 

        Parameters:
        -----------
        1) first_delta_bin:  (str) The first delta score bin to compare
        2) second_delta_bin: (str) The second delta score bin to compare

        Returns:
        ++++++++
        1) (str) The best/highest delta score bin
        """

        max_d_score1 = float(first_delta_bin.strip().split("-")[1])
        max_d_score2 = float(second_delta_bin.strip().split("-")[1])

        if float(max_d_score1) >= float(max_d_score2):
            return(first_delta_bin)
        else:
            return(second_delta_bin)

    ## error tracking
    cur_ref_allele_n = 0
    cur_empty_position_info = 0
    considered_position = set()
    added_var_pos = set()
    per_pos_max_delta_bin = defaultdict(lambda: delta_score_bins[0])

    ## Parse gnomad and get a dict of gnomAD varaints that pass all filters
    (gnomad_var_dict,
    non_matching_ref_alleles,
    non_matching_symbols) = parse_gnomad_by_region(region_chrom         = region_chrom,
                                                   region_start         = region_start,
                                                   region_end           = region_end,
                                                   region_strand        = region_strand,
                                                   contains_chr         = contains_chr,
                                                   vcf                  = vcf,
                                                   vep_index            = vep_index, 
                                                   spliceai_region_dict = spliceai_region_dict,
                                                   by_pos_cov_dict      = by_pos_cov_dict,
                                                   coverage_cutoff      = coverage_cutoff,
                                                   alt_symbol_dict      = alt_symbol_dict,
                                                   error_log_file       = error_log_file,
                                                   delta_score_bins     = delta_score_bins,
                                                   SCORE_TYPE           = SCORE_TYPE)


    ## Iterate through each position in the region and update O and E counts
    for region_pos in spliceai_region_dict.keys():

        ## If the region position does not overlap the query region, skip it
        if region_pos < region_start or region_pos > region_end:
            continue

        pos_list = []
        ## Check if position has a gnomAD variant
        if region_pos in gnomad_var_dict:
            
            ## Update the position list with the gnomAD variant info
            pos_list = gnomad_var_dict[region_pos]

        else: ## If no variant, get info for the variant free nucleotide
            
            ## Get ref allele
            ref_allele = list(spliceai_region_dict[region_pos].keys()).pop() if region_pos in spliceai_region_dict and len(spliceai_region_dict[region_pos].keys()) == 1 else "None"

            ## Get strand corrected alleles
            strand_corrected_ref = correct_allele_by_strand(region_strand, ref_allele)

            ## Alt allele is same as ref (No ref to alt change)
            strand_corrected_alt = strand_corrected_ref

            ## Skip any SpliceAI predictions with a ref allele of N
            if strand_corrected_ref == "N":
                error = ("\n\n!!ERROR!! Position with Ref allele 'N'"
                        " Position will be skpped"
                        " Position = {}:{}\n").format(region_chrom, region_pos)
                output_log(error,error_log_file)
                ref_allele_n += 1
                cur_ref_allele_n += 1
                continue

            ## skip sites with no info
            ### This happends when multiple ref allels occur but no ref allele matches the fasta ref allele
            if len(spliceai_region_dict[region_pos]) == 0:
                error = ("\n\n!!ERROR!! Missing Spliceai Info"
                        " at Position = {}:{}\n"
                        " Skipping Site").format(region_chrom, region_pos)
                output_log(error,error_log_file)
                empty_position_info += 1
                cur_empty_position_info += 1
                continue

            ## Get the max/sum of delta score among all ref to alt combinations for the current position 
            max_sum_delta_score = (max([get_max_sum_delta_score(spliceai_region_dict[region_pos][ref_allele][alt_allele]) for alt_allele in spliceai_region_dict[region_pos][ref_allele].keys()])
                                   if SCORE_TYPE == "sum" 
                                   else max([get_max_delta_score(spliceai_region_dict[region_pos][ref_allele][alt_allele]) for alt_allele in spliceai_region_dict[region_pos][ref_allele].keys()])
                                   if SCORE_TYPE == "max"
                                   else None)

            if max_sum_delta_score is None:
                raise ValueError("!!ERROR!! Unrecognized value for the SpliceAI Score Type == '{}'".format(SCORE_TYPE))

            delta_score_bin = get_delta_score_bin(max_sum_delta_score, delta_score_bins)

            zerotons = 1.0
            zerotons_plus_singletons = 1.0
            singletons = 0.0
            non_zerotons = 0.0
            non_zerotons_plus_singletons = 0.0 

            ## Update the position list
            pos_list.append({"zerotons": zerotons,
                             "zerotons_plus_singletons": zerotons_plus_singletons,
                             "singletons": singletons,
                             "non_zerotons": non_zerotons,
                             "non_zerotons_plus_singletons":non_zerotons_plus_singletons,
                             "strand_corrected_ref": strand_corrected_ref,
                             "strand_corrected_alt": strand_corrected_alt,
                             "max_sum_delta_score": max_sum_delta_score,
                             "delta_score_bin": delta_score_bin})


        ## Iterate through each dict in the position list
        for pos_dict in pos_list:
            
            ## Identify values
            score_bin = pos_dict["delta_score_bin"]
            corrected_ref = pos_dict["strand_corrected_ref"]
            corrected_alt = pos_dict["strand_corrected_alt"]
            nz = float(pos_dict["non_zerotons"])
            nzps = float(pos_dict["non_zerotons_plus_singletons"])

            ## Get zeroton and zeroton plus singleton mutation rate for the current delta score bin and ref allele
            zeroton_mutation_rate = float(mut_rate_dict[score_bin]["zeroton"][corrected_ref])
            zeroton_plus_singleton_mutation_rate = float(mut_rate_dict[score_bin]["zeroton_plus_singleton"][corrected_ref]) 
            
            ## add ovserved and expected scores
            ### Observed scores
            t_id_dict[transcript_id]["zeroton_observed_var_count"] += nz
            t_id_dict[transcript_id]["zero_and_singleton_observed_var_count"] += nzps

            t_id_dict[transcript_id]["{}_zeroton_observed".format(score_bin)] += nz
            t_id_dict[transcript_id]["{}_{}_zeroton_observed".format(score_bin, corrected_ref)] += nz
            t_id_dict[transcript_id]["{}_zero_and_singleton_observed".format(score_bin)] += nzps 
            t_id_dict[transcript_id]["{}_{}_zero_and_singleton_observed".format(score_bin, corrected_ref)] += nzps 

            ## Expected Scores
            t_id_dict[transcript_id]["zeroton_expectation_sum"] += zeroton_mutation_rate
            t_id_dict[transcript_id]["zero_and_singleton_expectation_sum"] += zeroton_plus_singleton_mutation_rate

            t_id_dict[transcript_id]["{}_zeroton_expected".format(score_bin)] += zeroton_mutation_rate 
            t_id_dict[transcript_id]["{}_{}_zeroton_expected".format(score_bin, corrected_ref)] += zeroton_mutation_rate 
            t_id_dict[transcript_id]["{}_zero_and_singleton_expected".format(score_bin)] += zeroton_plus_singleton_mutation_rate 
            t_id_dict[transcript_id]["{}_{}_zero_and_singleton_expected".format(score_bin, corrected_ref)] += zeroton_plus_singleton_mutation_rate

            ## Get the best delta score bin for the current position
            per_pos_max_delta_bin[region_pos] = get_best_delta_bin(per_pos_max_delta_bin[region_pos], score_bin) 

        
        ## Update considered positions
        considered_position.add(region_pos)
        
    ## Update the transcript dict with error info
    t_id_dict[transcript_id]["non_matching_ref_alelles"] += non_matching_ref_alleles
    t_id_dict[transcript_id]["non_matching_symbol"] += non_matching_symbols
    t_id_dict[transcript_id]["ref_alleles_as_N"] += cur_ref_allele_n
    t_id_dict[transcript_id]["empty_position"] += cur_empty_position_info 
    t_id_dict[transcript_id]["positions_considered"] += len(considered_position)


    ## Get the count of positions with each delta score bin
    delta_bin_counts = defaultdict(int)
    for pos,ds_bin in per_pos_max_delta_bin.items():
        delta_bin_counts[ds_bin] += 1

    ## Add delta score bin counts 
    for ds_bin, count in delta_bin_counts.items():
        t_id_dict[transcript_id]["{}_count".format(ds_bin)] += count

    return(non_matching_ref_alleles, 
           non_matching_symbols, 
           ref_allele_n, 
           empty_position_info)


def prepare_region(transcript_region_dict, 
                   tx_key, 
                   region_type,
                   window_size,
                   step_size): 
    """
    prepare_region
    ==============
    Method to create and prepare a gene/region for O and E processing. This method will create a dictionary used to 
     track O and E counts for a gene/region. 

    If the region type is set as 'gene', only one region will be created based on the gene's start and end position. 

    If the region type is set to "region", then mutliple regions will be generated for the given gene using the user 
     defined window size and step size. 

    Parameters:
    -----------
    1) transcript_region_dict: (dict) A dictionary for a specific gene/transcript generated by the 'get_transcript_regions' method. 
    2) tx_key:                  (str) The current gene/transcript name/key 
    3) region_type:             (str) The --region-type parameter. 'gene' or 'region'
    4) window_size:             (int) The intragenic window size to create region with. (Only used if region_type is set to "region")
    5) step_size:               (int) The step size to move the window along the gene. (Only used if the region_type is set to "region")

    Returns:
    ++++++++
    (dict) A dictionary where each top level key represents a region for the current gene to calculate O and E scores for. 
            lower level keys are used for data tracking for each region in the dict.
    """
        
    ## Create a region per gene dict
    regional_transcript_dict = defaultdict(lambda: defaultdict(str))

    if region_type == "region":
        regions = sliding_window_regions(start        = transcript_region_dict["gene_start"],
                                         end          = transcript_region_dict["gene_end"],
                                         window_size  = window_size, 
                                         step_size    = step_size) 

    elif region_type == "gene":
        regions = [[transcript_region_dict["gene_start"], transcript_region_dict["gene_end"]]]
    
    else:
        raise ValueError("!!ERROR!! Unrecognized input value for the --region-type argument '{}'".format(args.region_type))

    ## Add a new transcript = region key with region start and region end to the dict
    for region_start, region_end in regions:
        
        ## Add region info
        region_key = "{}:{}-{}".format(tx_key, region_start, region_end)
        regional_transcript_dict[region_key] = copy.deepcopy(transcript_region_dict)
        regional_transcript_dict[region_key]["start"] = region_start
        regional_transcript_dict[region_key]["end"] = region_end

        ## Add error tracking per transcript
        regional_transcript_dict[region_key]["non_matching_ref_alelles"] = 0.0
        regional_transcript_dict[region_key]["non_matching_symbol"] = 0.0
        regional_transcript_dict[region_key]["multiple_ref_allele_pos"] = 0.0
        regional_transcript_dict[region_key]["ref_alleles_as_N"] = 0.0
        regional_transcript_dict[region_key]["empty_position"] = 0.0
        regional_transcript_dict[region_key]["PAR_pos"] = 0.0
        regional_transcript_dict[region_key]["segdup_selfchain_repeat_pos"] = 0.0
        regional_transcript_dict[region_key]["positions_considered"] = 0.0
        regional_transcript_dict[region_key]["total_gene_positions"] = regional_transcript_dict[region_key]["gene_end"] - regional_transcript_dict[region_key]["gene_start"] + 1
        regional_transcript_dict[region_key]["total_region_positions"] = region_end - region_start + 1

    ## return new dict
    return(regional_transcript_dict)


def add_region_score_tracking(regional_transcript_dict, tx_key, score_bin_list):
    """
    add_region_score_tracking
    =========================
    Method to update the regional transcript dict with keys for tracking Observed and Expected counts.
     This method will add the updated keys to only the region designated by the tx_key. All other regions
     will be left alone. Keys are based on the SpliceAI score bins designated in the config file.

    Parameters:
    -----------
    1) regional_transcript_dict: (dict) The dictionary that contains the current regions to calculate O and E counts for. 
    2) tx_key:                    (str) The key in the dictionary for which region to add default count info for
    3) score_bin_list:           (list) A list of SpliceAI score bins to add to the region count tracker. (These bins are from the config file)

    Returns:
    ++++++++
    None: The dictionary is passed by reference, and therefore will be updated without need to return it.
    """

    ## Add tracking info
    for label in score_bin_list:

        ## Add label specific observed and expected score counts
        ### set observed score to 1 as default 
        ### Get the average mutation rate for the current bin and set as default for expectation   
        regional_transcript_dict[tx_key]["{}_zeroton_observed".format(label)] = 1.0
        regional_transcript_dict[tx_key]["{}_zero_and_singleton_observed".format(label)] = 1.0
        regional_transcript_dict[tx_key]["{}_zeroton_expected".format(label)] = 1.0
        regional_transcript_dict[tx_key]["{}_zero_and_singleton_expected".format(label)] = 1.0

        ## Delta score bin count
        regional_transcript_dict[tx_key]["{}_count".format(label)] = 0.0

        ## For each label, get a count per ref allele
        for ref_allele in ["A","C","G","T"]:
            
            ## Set observed to 1
            ## Set expected to the mutation rate for that ref allele
            regional_transcript_dict[tx_key]["{}_{}_zeroton_observed".format(label,ref_allele)] = 1.0
            regional_transcript_dict[tx_key]["{}_{}_zero_and_singleton_observed".format(label,ref_allele)] = 1.0
            regional_transcript_dict[tx_key]["{}_{}_zeroton_expected".format(label,ref_allele)] = 1.0
            regional_transcript_dict[tx_key]["{}_{}_zero_and_singleton_expected".format(label,ref_allele)] = 1.0  

    ## Add general observed and expected score counts
    regional_transcript_dict[tx_key]["zeroton_observed_var_count"] = 0.0
    regional_transcript_dict[tx_key]["zero_and_singleton_observed_var_count"] = 0.0
    regional_transcript_dict[tx_key]["zeroton_expectation_sum"] = 0.0
    regional_transcript_dict[tx_key]["zero_and_singleton_expectation_sum"] = 0.0


def create_output_file(out_file, score_bin_list):
    """
    create_output_file
    ===================
    Mehtod to create the regional score output file and the log output and error files. 
    """

    ## Create log files
    date_string = datetime.datetime.now().strftime("%m-%d-%Y")
    time_string = datetime.datetime.now().strftime("%H:%M:%S")
    output_log("##Start Date: {}\n##Start Time: {}\n#ERRORS:\n".format(date_string,time_string), 
               error_log_file,
               new_file = True)
              
    output_log("##Start Date: {}\n##Start Time: {}\n#STDOUT:\n".format(date_string, time_string), 
               out_log_file,
               new_file = True)

    ## Create file with header
    with open(out_file, "w") as out: 
        out_header = ["#chrom",
                      "region_start",
                      "region_end",
                      "txStart",
                      "txEnd",
                      "strand",
                      "gene_id",
                      "gene_symbol",
                      "transcript_id",
                      "max_exon_number",
                      "cds_exon_count",
                      "zeroton_observed_var_count",
                      "zeroton_expectation_sum",
                      "zero_plus_singleton_observed_var_count",
                      "zero_plus_singleton_expectation_sum"]

        out_header.extend(["{}_zeroton_observed".format(label) for label in score_bin_list])
        out_header.extend(["{}_zeroton_expected".format(label) for label in score_bin_list])
        out_header.extend(["{}_zero_and_singleton_observed".format(label) for label in score_bin_list])
        out_header.extend(["{}_zero_and_singleton_expected".format(label) for label in score_bin_list])
        out_header.extend(["{}_count".format(label) for label in score_bin_list])

        for d_bin in score_bin_list:
            for ref in ["A","C","G","T"]:
                
                out_header.append("{}_{}_zeroton_observed".format(d_bin,ref))
                out_header.append("{}_{}_zeroton_expected".format(d_bin,ref))
                out_header.append("{}_{}_zero_and_singleton_observed".format(d_bin,ref))
                out_header.append("{}_{}_zero_and_singleton_expected".format(d_bin,ref))

        out_header.append("non_matching_ref_alelles")
        out_header.append("non_matching_symbol")
        out_header.append("multiple_ref_allele_pos")
        out_header.append("ref_alleles_as_N") 
        out_header.append("empty_position")
        out_header.append("PAR_pos")
        out_header.append("segdup_selfchain_repeat_pos")
        out_header.append("positions_considered")
        out_header.append("total_gene_positions")
        out_header.append("total_region_positions")
        out.write("\t".join(out_header) + "\n")


def write_scored_region(out_file, regional_score_dict, transcript_key, score_bin_list):
    """
    write_scored_region
    ===================
    Method to write a scored region to the output file. This method takes the regional dict with the scored regions,
     removes the scored region from the dict, and writes it to the output file.

    NOTE: This method will remove the designated scored region from the regional score dict. 

    Parameters:
    -----------
    1) out_file: (str) The name (and path) of the output file to write to
    2) regional_score_dict: (dict) The dictionary with the scored region to write to the output file
    3) transcript_key:       (str) The unique regional transcript key to use to identify which region to write to the output file
    4) score_bin_list:      (list) A list of SpliceAI score bins from the config file

    Returns:
    ++++++++
    None: This function does not return anything. It only writes the specific scored region to the output file. The scored region 
           that is written to the output file will be permanently removed from the regional score dict. 
    """
            
    ## Write region specific information to file
    with open(out_file, "a") as out: 

        ## Get the value for the current region and remove that region from the dict
        value = regional_score_dict.pop(transcript_key)

        output_list = [str(value["chrom"]), 
                       str(value["start"]), 
                       str(value["end"]), 
                       str(value["gene_start"]),
                       str(value["gene_end"]),
                       str(value["strand"]), 
                       str(value["gene_id"]), 
                       str(value["gene_name"]), 
                       str(value["transcript_id"]),
                       str(value["max_exon_number"]),
                       str(value["cds_exon_count"]),
                       str(value["zeroton_observed_var_count"]), 
                       str(value["zeroton_expectation_sum"]),
                       str(value["zero_and_singleton_observed_var_count"]),
                       str(value["zero_and_singleton_expectation_sum"])
                       ]    

        output_list.extend([str(value["{}_zeroton_observed".format(label)]) for label in score_bin_list])
        output_list.extend([str(value["{}_zeroton_expected".format(label)]) for label in score_bin_list])
        output_list.extend([str(value["{}_zero_and_singleton_observed".format(label)]) for label in score_bin_list])
        output_list.extend([str(value["{}_zero_and_singleton_expected".format(label)]) for label in score_bin_list])
        output_list.extend([str(value["{}_count".format(label)]) for label in score_bin_list])

        for d_bin in score_bin_list:
            for ref in ["A","C","G","T"]:
                
                output_list.append(str(value["{}_{}_zeroton_observed".format(d_bin,ref)]))
                output_list.append(str(value["{}_{}_zeroton_expected".format(d_bin,ref)]))
                output_list.append(str(value["{}_{}_zero_and_singleton_observed".format(d_bin,ref)]))
                output_list.append(str(value["{}_{}_zero_and_singleton_expected".format(d_bin,ref)]))


        output_list.append(str(value["non_matching_ref_alelles"]))
        output_list.append(str(value["non_matching_symbol"]))
        output_list.append(str(value["multiple_ref_allele_pos"]))
        output_list.append(str(value["ref_alleles_as_N"]))
        output_list.append(str(value["empty_position"]))
        output_list.append(str(value["PAR_pos"]))
        output_list.append(str(value["segdup_selfchain_repeat_pos"]))
        output_list.append(str(value["positions_considered"]))
        output_list.append(str(value["total_gene_positions"]))
        output_list.append(str(value["total_region_positions"]))

        out.write("\t".join(output_list) + "\n")
        output_list = []


#---------------------------------------------------------------------------------------------------------------------------------
## Main
#---------------------------------------------------------------------------------------------------------------------------------

def o_e_counts(parser, args):
    """
    o_e_counts 
    ==========
    Main method used to control the calculation of the Observed and Expected counts for the ConSplice Model.
    """

    print("\n\n\t********************************************")
    print("\t* ConSplice - Observed and Expected Counts *")
    print("\t********************************************\n\n")

    start_time = time.time()

    if args.chrom == None:
        args.chrom = ["All"]


    ## Check args
    if args.region_type == "region":

        if args.window_size is None:
            print("\n!!ERROR!! The --window-size argument is required when --region-type is set to 'region'.\n")
            sys.exit(1)

        if args.step_size is None:
            print("\n!!ERROR!! The --step-size argument is required when --region-type is set to 'region'.\n")
            sys.exit(1)


        if int(args.window_size) < 1:
            print("\n!!ERROR!! The window size needs to be a positive value. Please enter a valid window size\n")
            sys.exit(1)

        if int(args.step_size) < 1:
            print("\n!!ERROR!! The step size needs to be a positive value. Please enter a valid step size\n")
            sys.exit(1)

    elif args.region_type == "gene":
        
        if args.window_size is not None:
            print("\n**WARNING** The --window-size argument was set when the --region-type was set to 'gene'. The --window-size will no be used\n")

        if args.step_size is not None:
            print("\n**WARNING** The --step-size argument was set when the --region-type was set to 'gene'. The --step-size will no be used\n")


    ## Check coverage cutoff value
    if args.cc > 1.0 or args.cc < 0.0: 
        print("\n!!ERROR!! The --cc (Coverage Cutoff) is set at {}. It needs to be between 0.0 and 1.0. Please correct the problem and try again".format(args.cc))
        sys.exit(1)

    ## Check number of CPUs
    ncpus = args.n_cpu if args.n_cpu > 0 else 1
    if ncpus > multiprocessing.cpu_count(): 
        print("\n**Warning** Too many CPUs designated. Number of CPUs will be set to {}".format( multiprocessing.cpu_count() if multiprocessing.cpu_count() > 0 else 1)) 
        ncpus = multiprocessing.cpu_count() if multiprocessing.cpu_count()  > 0 else 1

    print(("\nInput Arguments:"
           "\n================"
           "\n - config-path:          {}"
           "\n - gnomad-vcf:           {}"
           "\n - coverage:             {}"
           "\n - gtf-file:             {}"
           "\n - spliceai-vcf:         {}"
           "\n - alt-gene-symbol:      {}"
           "\n - fasta:                {}"
           "\n - seg-dups              {}"
           "\n - self-chains           {}" 
           "\n - out-file:             {}"
           "\n - substitution-matrix:  {}"
           "\n - region-type:          {}"
           "\n - window-size:          {} {}"
           "\n - step-size:            {} {}" 
           "\n - chrom:                {}"
           "\n - cc:                   {}"
           "\n - coverage-label:       {}"
           "\n - n-cpu:                {}"
           "\n - repeat-score-cutoff   {}"
           "\n - spliceai-score-type:  {}"
           ).format(args.config_path,
                    args.gnomad_vcf, 
                    args.coverage, 
                    args.gtf_file, 
                    args.spliceai_vcf, 
                    args.alt_gene_symbol,
                    args.fasta, 
                    args.seg_dups,
                    args.self_chains,
                    args.out_file, 
                    args.substitution_matrix,
                    args.region_type, 
                    args.window_size, "(Not used)" if args.region_type == "gene" else "",
                    args.step_size, "(Not used)" if args.region_type == "gene" else "",
                    ", ".join(args.chrom), 
                    args.cc, 
                    args.coverage_label,
                    args.n_cpu,
                    args.repeat_score_cutoff,
                    args.spliceai_score_type,
                    )
    )


    print("\n\nSetting up ConSplice Configs")
    print("============================")

    load_config_file(config_path = args.config_path, 
                     region_type = args.region_type, 
                     chrom_list = args.chrom, 
                     spliceai_score_type = args.spliceai_score_type,
                     window_size = args.window_size if args.window_size is not None else "GENE",
                     step_size = args.step_size if args.step_size is not None else "GENE")

    
    print("\n\nGetting gene and transcript region info")
    print("=======================================")

    print("\n\tGathering all transcript regions")
    transcript_region_dict, gene_strand_dict = get_transcript_regions(args.gtf_file, biotype="protein_coding")

    print("\t\tCollected {} transcripts for processing".format(len(transcript_region_dict)))

    print("\n\tCreating an map of alternative gene symbols")
    ## Get a dictionary of alternative gene symbols
    alt_symbol_dict = get_alternative_gene_symbols(args.alt_gene_symbol)

    print("\n\tUpdating gene symbol mapping with strand info")
    ## Update gene strand dict with alt symbols 
    gene_strand_dict = expand_gene_dict(alt_symbol_dict, gene_strand_dict)


    print("\nCreating Interlap Objects for exclude regions") 
    print("=============================================") 

    ## Create an interlap object for seg dups and for self chains. 
    print("\n\tCreating Interlap Object for Segmental Duplications")
    segdup_interlap_dict = create_interlap_from_ggd_pkg(file_path = args.seg_dups, 
                                                        zero_based = True, 
                                                        set_cutoff = True, 
                                                        cutoff_col = "fracMatch", 
                                                        cutoff_score = args.repeat_score_cutoff )

    print("\n\tCreating Interlap Object for High Identity Self Chains")
    selfchain_interlap_dict = create_interlap_from_ggd_pkg(file_path = args.self_chains, 
                                                           zero_based = True, 
                                                           set_cutoff = True, 
                                                           cutoff_col = "normalized_score", 
                                                           cutoff_score = (args.repeat_score_cutoff * 100), 
                                                           chrom_col = "target_chrom", 
                                                           start_col = "target_start", 
                                                           end_col = "target_end")

    ## Get the index of the row that represents the header
    print("\n\nGetting the substitution rates from the substitution matrix")
    print("===========================================================")
    header_index = 0
    try:
        with io.open(args.substitution_matrix, "rt", encoding = "utf-8") as mrt:
            for i,line in enumerate(mrt):
                if line[0] == "#": 
                    header_index = i
                elif line[0] != "#":
                    break
    except IOError as e:
        print("\n\n!!ERROR!! Unable to read the substitution matrix file")
        print(str(e))
        sys.exit(1)

    print("\n\tUsing the SpliceAI aware subtitution matrix: '{}'".format(args.substitution_matrix))
    ## Load the matrix into a pandas df
    mr_table = pd.read_csv(args.substitution_matrix, sep="\t", index_col=False, header = header_index)
    mr_table = mr_table.rename(columns = {"#ref":"ref"})


    print("\n\tGetting substitution rates")
    ## Get the substitution rate for each SpliceAI bin and ref allele combo
    by_ref_mutation_rate_dict = defaultdict(lambda: defaultdict(float))
    mr_table.groupby(["delta_score_bin"]).apply(get_normalized_matrix_per_label, mutation_rate_dict = by_ref_mutation_rate_dict)

    print("\n\tChecking that SpliceAI score bins from the config file match the substitution rate bins")
    
    assert sorted(delta_score_bins) == sorted(list(by_ref_mutation_rate_dict.keys())), ("\n!!ERROR!! The Substitution rate bins don't match the SpliceAI Score Bins in the config file"
                                                                                        "\n\t- Substitution rate bins: {}"
                                                                                        "\n\t- SpliceAI Score bins:    {}\n").format( ", ".join(sorted(list(by_ref_mutation_rate_dict.keys()))), 
                                                                                                                                      ", ".join(sorted(delta_score_bins)))

    print("\nSubstitution rates by SpliceAI bin and reference allele")
    print("-------------------------------------------------------")
    for label,value in by_ref_mutation_rate_dict.items():
        print("SpliceAI score bin: {}".format(label))
        for zeroton,value2 in value.items():
            print("\tRates for {} model".format(zeroton))
            for ref,mr in value2.items():
                print("\t\t",ref + ":", mr)
            print("\n")
        print("\n")
    print("-------------------------------------------------")

    print("\n\nPreparing VCF parsing")
    print("=====================")

    print("\n\tLoading the gnomAD vcf file")
    ## load the gnomAD vcf file using cyvcf2
    gnomad_vcf = VCF(args.gnomad_vcf,threads=ncpus)

    ## Check if the gnmoad vcf file has a chr prefix or not
    contains_chr = True if any("chr" in x for  x in gnomad_vcf.seqnames) else False

    ## Get the VEP info field and create a list for indexing 
    gnomad_vep_index = gnomad_vcf["vep"]["Description"].split(":")[-1].strip().replace('"',"").split("|") 


    print("\n\tLoading the gnomAD coverage file")
    try:
        ## Load tabixed coverage file
        cov_file = pysam.TabixFile(args.coverage)
       
       ## Get the header index from the coverage file header
        cov_header = cov_file.header[0].strip().split("\t")
        cov_header_dict = {x:i for i,x in enumerate(cov_header)}

    except IOError as e:
        print(str(e))
        print("\n!!ERROR!! Please provide a bgzipped and tabixed coverage file file\n")
        sys.exit(1)
    

    print("\n\tLoading SpliceAI vcf file")
    ## load the SpliceAI vcf file using cyvcf2
    spliceai_vcf = VCF(args.spliceai_vcf, threads=ncpus)

    ## Get the SpliceAI info and create an index
    spliceai_info_index = spliceai_vcf["SpliceAI"]["Description"].strip().split("Format:")[1].strip().replace("\"","").split("|")

    print("\n\nPreapting transcript variant tracking")
    print("=====================================")

    delta_score_bin_list = list(by_ref_mutation_rate_dict.keys())

    create_output_file(out_file = args.out_file, score_bin_list = delta_score_bin_list)

    if args.region_type == "region":

        print("\n\tSplitting genes into regions using a sliding window and preparing region tracking")
        print("\n\tWindow Size = {}".format(args.window_size))
        print("\n\tStep Size = {}".format(args.step_size))

    elif args.region_type == "gene":

        print("\n\tUsing gene positions as constraint regions")


    prev_chrom = 1

    ## Track number of non matching sites
    non_matching_ref_alleles = 0
    non_matching_symbols = 0
    pos_with_multiple_refs = 0
    ref_allele_n = 0
    empty_position_info = 0
    seg_dup_or_self_chain = 0


    ## Iterate through the transcripts
    print("\n\nParsing transcripts and calculating O and E counts")
    print("==================================================\n")
    for i,tx_key in enumerate(transcript_region_dict.keys()):

        if i % 500 == 0:
            print("\t- Processed {} transcripts".format(i))
            output_log("\n\t- Processed {} transcripts".format(i),out_log_file)

        ## Create a dictionary to track the regions for the current gene
        regional_score_dict = prepare_region(transcript_region_dict = transcript_region_dict[tx_key], 
                                                         tx_key = tx_key, 
                                                         region_type = args.region_type,
                                                         window_size = args.window_size,
                                                         step_size = args.step_size)

        ## iterate over each region and get o and e scores
        dict_keys = list(regional_score_dict.keys())
        for transcript_key in dict_keys:

            ## Add score tracking info
            add_region_score_tracking(regional_transcript_dict = regional_score_dict, 
                                      tx_key = transcript_key,
                                      score_bin_list = delta_score_bin_list)


            ## Get region info
            chrom = regional_score_dict[transcript_key]["chrom"]
            start = regional_score_dict[transcript_key]["start"]
            end = regional_score_dict[transcript_key]["end"]
            strand = regional_score_dict[transcript_key]["strand"]
            gene_name = regional_score_dict[transcript_key]["gene_name"] 

            ## Skip any non selected chromosomes 
            if args.chrom != ["All"] and chrom not in args.chrom:
                continue

            ## reload vcf
            if prev_chrom != chrom:
                spliceai_vcf = VCF(args.spliceai_vcf, threads=ncpus)
                gnomad_vcf = VCF(args.gnomad_vcf, threads=ncpus)
                prev_chrom = chrom

            ## Check for repeats in the current region
            nonrepeat_regions, repat_count = check_repeats(chrom, 
                                                           start, 
                                                           end, 
                                                           segdup_interlap_dict[chrom], 
                                                           selfchain_interlap_dict[chrom])

            seg_dup_or_self_chain += repat_count
            regional_score_dict[transcript_key]["segdup_selfchain_repeat_pos"] += repat_count  

            ## Check the PAR region for the current region if on X chromosome
            if chrom == "X":

                new_regions = []
                ## Iterate over the non repeat regions
                for nr_region in nonrepeat_regions:

                    ## if the non repeat region is empty, meaning the entire region was overlapped by a repeat, skip this region
                    if nr_region.size == 0:
                        continue 

                    ## Get all non par regions as a numpy array
                    non_par, par_count = check_xchrom_par(chrom, nr_region[0], nr_region[1], PAR_PR)

                    ## Update PAR bases for the current region
                    regional_score_dict[transcript_key]["PAR_pos"] += par_count 

                    ## If the non par array is not empty add it to the new/filtered regions
                    if non_par.size > 0:
                        new_regions.append(non_par)

                ## Update that non repeat regions array with the filtered non PAR regions
                nonrepeat_regions = np.concatenate(new_regions) if len(new_regions) > 0 else np.array([[]])


            ## Get SpliceAI region Info
            spliceai_region_list, spliceai_position_dict, multi_ref_pos = spliceai_scores_by_region(region_list = nonrepeat_regions,
                                                                                                    spliceai_vcf = spliceai_vcf,
                                                                                                    spliceai_info_index = spliceai_info_index,
                                                                                                    chrom = chrom,
                                                                                                    gene_name = gene_name,
                                                                                                    strand = strand,
                                                                                                    alt_symbol_dict = alt_symbol_dict,
                                                                                                    gene_strand_dict = gene_strand_dict,
                                                                                                    check_symbol = True,
                                                                                                    fasta_file = args.fasta,
                                                                                                    error_log_file = error_log_file)


            ## out log
            stdout_string = ("\nBefore Filtering Region: {}:{}-{}"
                             "\nSpliceAI Filtered Region(s): {}").format(chrom,start,end,spliceai_region_list)
            output_log(stdout_string,out_log_file)

            ## the number of position with multiple ref alleles 
            regional_score_dict[transcript_key]["multiple_ref_allele_pos"] = len(multi_ref_pos)

            ## Iterate over each relateive spliceai continous region
            for relative_spliceai_region in spliceai_region_list:

                ## Skip regions where no spliceai info was found
                if relative_spliceai_region["start"] == float("inf") and relative_spliceai_region["end"] == 0:
                    continue

                ## Get a dict of coverage where key = position, value = coverage value
                ### {POS:COVERAGE,POS:COVERAGE,POS:COVERAGE,...}
                by_position_coverage_dict = gnomad_coverage_by_region(cov_header_dict = cov_header_dict, 
                                                                      coverage_label = args.coverage_label, 
                                                                      cov_file = cov_file, 
                                                                      region = relative_spliceai_region) 


                ## Update mutation frequencies from gnomAD variants
                (non_matching_ref_alleles, 
                 non_matching_symbols, 
                 ref_allele_n, 
                 empty_position_info) = get_by_region_o_e_counts(gnomad_vcf,
                                                                    contains_chr,
                                                                    relative_spliceai_region["chrom"],
                                                                    relative_spliceai_region["start"],
                                                                    relative_spliceai_region["end"],
                                                                    strand,
                                                                    by_position_coverage_dict,
                                                                    args.cc,
                                                                    gnomad_vep_index,
                                                                    spliceai_position_dict,
                                                                    alt_symbol_dict,
                                                                    regional_score_dict,
                                                                    transcript_key,
                                                                    by_ref_mutation_rate_dict,
                                                                    non_matching_ref_alleles,
                                                                    non_matching_symbols,
                                                                    ref_allele_n,
                                                                    empty_position_info)
            

            ## Write the scored region to the output file
            write_scored_region(out_file = args.out_file,
                                regional_score_dict = regional_score_dict,
                                transcript_key = transcript_key,
                                score_bin_list = delta_score_bin_list)

    end_time = time.time()

    print("\t- Processed {} transcripts".format(i + 1))
    output_log("\n\t- Processed {} transcripts".format(i + 1),out_log_file)

    ## Update log files
    output_log("\nElapsed Time in Seconds: {}".format(end_time - start_time) ,out_log_file)
    output_log("\n\nNumber of sites where the reference alleles between SpliceAI and gnomAD did not match: {}".format(non_matching_ref_alleles), error_log_file)
    output_log("Number of sites where the gene symbols did not match between SpliceAI and gnomAD: {}".format(non_matching_symbols), error_log_file)
    output_log("Number of positions where SpliceAI had mutliple reference alleles: {}".format(pos_with_multiple_refs), error_log_file)
    output_log("Number of positions where SpliceAI had a reference allele 'N': {}".format(ref_allele_n), error_log_file)
    output_log("Number of positions where SpliceAI had multiple ref alleles, all of which did not match the fasta ref alelle: {}".format(empty_position_info), error_log_file)
    output_log("Number of positions removed because of Seg Dup or Self Chain overlap: {}".format(seg_dup_or_self_chain), error_log_file)

    ## Print report 
    print("\n\n========================================================================")
    print("General STDOUT:")
    print("---------------")
    print("STDOUT LOG written to {}".format(out_log_file))
    print("\nERROR Report:")
    print("-------------")
    print(" Number of sites where the reference alleles between SpliceAI and gnomAD did not match: {}".format(non_matching_ref_alleles))
    print(" Number of sites where the gene symbols did not match between SpliceAI and gnomAD: {}".format(non_matching_symbols))
    print(" Number of positions where SpliceAI had mutliple reference alleles: {}".format(pos_with_multiple_refs))
    print(" Number of positions where SpliceAI had a reference allele 'N': {}".format(ref_allele_n))
    print(" Number of positions where SpliceAI had multiple ref alleles, all of which did not match the fasta ref alelle: {}".format(empty_position_info))
    print(" Number of positions removed because of Seg Dup or Self Chain overlap: {}".format(seg_dup_or_self_chain))
    print(" ERROR LOG written to {}".format(error_log_file))
    print("\nElapsed Time in Seconds: {}".format(end_time - start_time))
    print("========================================================================")


    print("\nOutput written to: '{}'".format(args.out_file))

    print("\nDONE")
