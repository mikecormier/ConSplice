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
from collections import defaultdict
from tqdm import tqdm
import json
import multiprocessing
import datetime
from pyfaidx import Fasta
import time
from .utils import (create_interlap_from_ggd_pkg,
                   correct_allele_by_strand, 
                   expand_gene_dict_with_alt_symbols, 
                   get_alternative_gene_symbols, 
                   get_delta_score_bin, 
                   get_max_delta_score, 
                   get_max_sum_delta_score,
                   output_log)


#---------------------------------------------------------------------------------------------------------------------------------
## Global Vars
#---------------------------------------------------------------------------------------------------------------------------------
delta_score_bins = ["0.0-0.01",
                    "0.01-0.05",
                    "0.05-0.1",
                    "0.1-0.25",
                    "0.25-0.5",
                    "0.5-1.0",
                    "1.0-1.75",
                    "1.75-4.0"]

## Log files
file_date = datetime.datetime.now().strftime("%m-%d-%Y")
error_log_file = file_date + "_sum_8_bin_o_and_e_scores.ERROR.log_count_matrix.ERROR.log"
out_log_file = file_date + "_sum_8_bin_o_and_e_scores.OUT.log"

#---------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
#---------------------------------------------------------------------------------------------------------------------------------

def add_gene_o_e(sub_p):

    p = sub_p.add_parser("gene-oe",
                         help = "Calculate Observed and Expected splicing variant counts for full gene regions",
                         formatter_class=argparse.RawDescriptionHelpFormatter,
                         description=("\n\t***********************\n"
                                      "\t* Gene O and E Scores *\n"
                                      "\t***********************\n\n"
                                      "\tCalculate Observed splicing variation from gnomAD and"
                                      "\n\tExpected splicing variation using a Substitution"
                                      "\n\tmatrix for full gene reiongs."
                        )
    )

    req = p.add_argument_group("Required Arguments")
    req.add_argument(
        "--gnomad-vcf",
        metavar = "gnomAD VCF/BCF File", 
        required=True, 
        help="(Required) The path to the gnomAD vcf/bcf file."
    )

    req.add_argument(
        "--coverage",
        metavar = "gnomAD Coverage file",
        required = True,
        help = "(Required) The gnomAD coverage file that goes with the gnomAD vcf file. The coverage files should be in a bed formated, bgzipped, and tabixed. (Use the ggd gnomad coverage recipe)"
    )

    req.add_argument(
        "--gtf-file",
        metavar = "GTF file",
        required=True,
        help="(Required) A gtf file to get gene genomic coordinance from to build the count matrix."
    )   

    req.add_argument(
        "--spliceai-vcf",
        metavar = "SpliceAI SNV prediction File", 
        required=True, 
        help="(Required) The path to SpliceAI delta score predictions for every snv vcf file."
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
        help="(Required) A fasta file to get reference alleles by position from"
    )   

    req.add_argument(
        "--seg-dups",
        metavar="Segmental Duplications File",
        required=True,
        help="(Required) The file path to the segmental duplications file. This file will be used to identify regions that should not be used for the calcualtion. (We suggest using the ggd segmental duplications data package)"
    )

    req.add_argument(
        "--self-chains",
        metavar="High Identity Self Chains",
        required=True,
        help="(Required) The file path to the high identity self chain repeat file. This file will be used to identify regions that should not be used for the calcualtion. (We suggest using the ggd self chain high identity data package)"
    )

    req.add_argument(
        "--out-file",
        metavar = "Output file",
        required = True,
        help = "(Required) The name of the output file to create."
    )

    req.add_argument(
        "-mrt",
        "--mutation-rate-table",
        metavar="SpliceAI Mutation Rate Table",
        required=True,
        help="Path and/or name of the mutation rate table created using SpliceAI predictions and gnomAD observed varaints. (Created from script: create_count_matrix.py"
    )

    p.add_argument(
        "--chrom",
        metavar = "Chromosomes",
        nargs = "*",
        choices = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","All"],
        help = "Which chromosomes to use to build the count matrix. For autosome only choose chromosomes 1-22, for X choose only X, All for all, etc., meaning all reference chromosomes"
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

    p.set_defaults(func=gene_o_e)


#---------------------------------------------------------------------------------------------------------------------------------
## Functions/Methods
#---------------------------------------------------------------------------------------------------------------------------------


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
    1) (dict) A dictionary of transcript regions. 1st key is the transcript id, 2nd layer keys include ["chrom","start","end","strand","gene_id","gene_name","max_exon_number", and "cds_exon_count"]  
    2) (dict) A dictionary of gene names (symbols) to gene ids 
    """
    ## Parse GTF file and create a pandas dataframe
    try:
        fh = gzip.open(gtf_file, "rt", encoding = "utf-8") if gtf_file.endswith(".gz") else io.open(gtf_file, "rt", encoding = "utf-8") 
        line_count = sum(1 for line in fh)
        fh.close()
    except IOError as e:
        print("!!Error!! Unable to read the gtf file: '{}'".format(gtf_file))
        print(str(e))
        sys.exit(1)

    ## Transcript dict
    transcript_region_dict = defaultdict(lambda: defaultdict(str))
    gene_name_dict = defaultdict(lambda: defaultdict(str))

    print("\n\tGetting {} transcript info from the gtf file".format(biotype))
    ## progress bar
    pbar = tqdm(total = line_count)
    pbar.set_description("Parsing GTF File")

    ## Iterate over gtf file
    try:
        fh = gzip.open(gtf_file, "rt", encoding = "utf-8") if gtf_file.endswith(".gz") else io.open(gtf_file, "rt", encoding = "utf-8") 
    except IOError as e:
        print("!!Error!! Unable to read the gtf file: '{}'".format(gtf_file))
        print(str(e))
        sys.exit(1)

    header = ["chrom","source","feature","start","end","score","strand","frame","attribute"]
    for line in fh:

        pbar.update(1)

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
                    transcript_region_dict[line_dict["transcript_id"]]["start"] = int(line_dict["start"])
                    transcript_region_dict[line_dict["transcript_id"]]["end"] = int(line_dict["end"])
                    transcript_region_dict[line_dict["transcript_id"]]["strand"] = line_dict["strand"]
                    transcript_region_dict[line_dict["transcript_id"]]["gene_id"] = line_dict["gene_id"]
                    transcript_region_dict[line_dict["transcript_id"]]["gene_name"] = line_dict["gene_name"]

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
    pbar.close()
                    
    return(transcript_region_dict, gene_name_dict)


def get_normalized_matrix_per_label(groupby_object, mutation_rate_dict):
    """
    get_normalized_matrix_per_label
    ===============================
    This method is used to get the marignal distribution for each reference allele at each near splice site bin (exp: D+5), 
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
    2) mutation_rate_dict: (dict)      A dictionary that contains the per nss position ref allele specific mutation rates. (Updated by this function)

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


def get_gnomad_vars_by_position(vcf,
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
    get_gnomad_vars_by_position
    ===========================
    This method is used to get the gnomAD variants between a genomic region (start and end position). Variants at any of the 
     positions in the region that pass filtering are used as an observed count for a gene at a specific bin. 

    It takes a cyvcf2 object and iterates over ecah varaint that is found between the position interval. Checks are done 
     for each varaint found to ensure the variant is appropriate. Checks include: 
        1) If the variant is an SNV. If the variant is not an SNV the variant will be skipped. 
        2) If the variant has a PASS filter. If the variant does not have a PASS filter the variant is skipped
        3) If the variant reference allele matches the reference allele anotated in the near splice site being compared
        4) If the variant has sufficient sequencing coverage. If not, the variant is skipped.
        5) If the variant is wihtin a transcript that the near splice site position is annotated in. If so, it is kept.
            if not, the variant is skipped. 
    
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
    1)  vcf:                      (cyvcf2 
                                  obj)    A cyvcf2 object for the current vcf/bcf file being queried 
    2)  contains_chr:             (bool)  whether or not the the sequence names in the vcf start with chr or not
    3)  region_chrom:             (str)   The chromsome to look in for the position interval 
    4)  region_start:             (int)   The start position of the interval
    5)  region_end:               (int)   The end position of the interval
    6)  region_strand:            (str)   A string representing the strand the gene is on. (+ or -)
    7)  by_pos_cov_dict:          (dict)  A dictionary of var coverage for the current region
    8)  coverage_cutoff:          (float) The coverage value at each position to use as a cutoff from the coverage file. (Between 0.0 and 1.0) 
    9)  vep_index:                (list)  A list that represents the vep csq index in the gnomad vcf 
    10) spliceai_region_dict:     (dict)  A dictionary of SpliceAI information for the current region 
    11) alt_symbol_dict:          (dict)  A dictionary of different alternative gene symbols 
    12) t_id_dict:                (dict)  A by transcript id dict that tracks the observed and expected scores for each transcript id
    13) transcript_id:            (str)   The transcript id for the current transcript being investigated
    14) mut_rate_dict             (dict)  A dictionary that contains the mutation rate information
    15) non_matching_ref_alleles: (int)   The current count of non matching ref alleles
    16) non_matching_symobls:     (int)   The current count of non matching symbols
    17) ref_allele_n:             (int)   The current count of reference alleles that are 'N'
    18) empty_position_info:      (int)   The current count of positions with null information 


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
    cur_non_matching_ref_alleles = 0
    cur_non_matching_symbols = 0
    cur_ref_allele_n = 0
    cur_empty_position_info = 0
    considered_position = set()
    added_var_pos = set()
    per_pos_max_delta_bin = defaultdict(lambda: delta_score_bins[0])


    ## Iterate over varaints that are with the same chromosome and position
    search = "chr{}:{}-{}".format(region_chrom,region_start,region_end) if contains_chr else "{}:{}-{}".format(region_chrom,region_start,region_end)  
    for var in vcf(search):

        ## Get variant info
        chrom = var.CHROM
        pos = var.POS
        ref = var.REF
        alt = var.ALT[0]    
        AC = var.INFO.get("AC")
        var_filter = var.FILTER if var.FILTER != None else "PASS"
        zerotons = 0
        zerotons_plus_singletons = 0
        singletons = 0
        non_zerotons = 0
        non_zerotons_plus_singletons = 0 

        ## Skip non SNVs
        ### We will only consider SNVs
        if len(ref) > 1 or len(alt) > 1:
            continue

        ## Only look at variants with a PASS filter
        if var_filter != "PASS":
            continue

        ## Check for matching reference alleles at the current position
        spliceai_ref_allele = list(spliceai_region_dict[pos].keys()).pop() if pos in spliceai_region_dict and len(spliceai_region_dict[pos].keys()) == 1 else "None" 
        if str(ref) != str(spliceai_ref_allele):
            ## Log error
            different_ref = list(spliceai_region_dict[pos].keys()).pop() if pos in spliceai_region_dict and len(spliceai_region_dict[pos].keys()) == 1 else "None"
            error = ("\n\n!!ERROR!! REF alleles don't match"
                     "\nPOSITION: {}"
                     "\ndiffering REF: {}"
                     "\nvariant info: {}"
                     "This position will be skipped").format(pos,different_ref, var)
            output_log(error,error_log_file)
            non_matching_ref_alleles += 1
            cur_non_matching_ref_alleles += 1

            added_var_pos.add(pos)
            continue

        ## Get the position specific coverage
        ## If position does not have enough sequence coverage, skip it 
        if float(by_pos_cov_dict[int(pos)]) < coverage_cutoff:
            continue

        ## Get a list of dictionaries for each vep annotation for the current var. 
        vep_annotations = [dict(zip(vep_index, x.strip().split("|"))) for x in var.INFO["vep"].strip().split(",")]

        ## Iterate over each vep annotation
        found = False
        var_symbols = []
        for csq in vep_annotations:

            ## Check for matching Symbols 
            if csq["SYMBOL"] in spliceai_region_dict[pos][ref][alt].keys() or csq["SYMBOL"] in { alt_symbol for symbol in spliceai_region_dict[pos][ref][alt].keys() if symbol in alt_symbol_dict for alt_symbol in alt_symbol_dict[symbol] }:
                found = True
                break
            else:
                var_symbols.append(csq["SYMBOL"])
                

        ## If the variant is unable to match with the spliceai, skip it
        if not found:
            error = ("\n\n**WARNING** NON matching gene symbols."
                     "\nPOSITION: {}"
                     "\n\tSite Gene: {}"
                     "\n\tVariant Gene: {}"
                     "\nSite will be skipped.").format(pos, ", ".join(spliceai_region_dict[pos][ref][alt].keys()), ", ".join(var_symbols))
            output_log(error,error_log_file)
            non_matching_symbols += 1
            cur_non_matching_symbols += 1

            added_var_pos.add(pos)
            continue

        ## Get strand corrected alleles
        strand_corrected_ref = correct_allele_by_strand(region_strand, ref)
        strand_corrected_alt = correct_allele_by_strand(region_strand, alt)

        ## Get the sum of delta scores for the current position, reference allele, and alternative allele
        max_sum_delta_score = get_max_sum_delta_score(spliceai_region_dict[pos][ref][alt])
        delta_score_bin = get_delta_score_bin(max_sum_delta_score, delta_score_bins)

        ## AC type
        zerotons = 0.0
        zerotons_plus_singletons = 1.0 if AC <= 1 else 0.0
        singletons = 1.0 if AC == 1 else 0.0
        non_zerotons = 1.0 if AC > 0.0 else 0.0
        non_zerotons_plus_singletons = 1.0 if AC > 1 else 0.0 

        ## Get zeroton and zeroton plus singleton mutation rate for the current delta score bin and ref allele
        zeroton_mutation_rate = float(mut_rate_dict[delta_score_bin]["zeroton"][strand_corrected_ref])
        zeroton_plus_singleton_mutation_rate = float(mut_rate_dict[delta_score_bin]["zeroton_plus_singleton"][strand_corrected_ref]) 


        ## add ovserved and expected scores
        ### Observed scores
        t_id_dict[transcript_id]["zeroton_observed_var_count"] += non_zerotons
        t_id_dict[transcript_id]["zero_and_singleton_observed_var_count"] += non_zerotons_plus_singletons

        t_id_dict[transcript_id]["{}_zeroton_observed".format(delta_score_bin)] += non_zerotons
        t_id_dict[transcript_id]["{}_{}_zeroton_observed".format(delta_score_bin,strand_corrected_ref)] += non_zerotons
        t_id_dict[transcript_id]["{}_zero_and_singleton_observed".format(delta_score_bin)] += non_zerotons_plus_singletons 
        t_id_dict[transcript_id]["{}_{}_zero_and_singleton_observed".format(delta_score_bin,strand_corrected_ref)] += non_zerotons_plus_singletons 

        ## Expected Scores
        t_id_dict[transcript_id]["zeroton_expectation_sum"] += float(zeroton_mutation_rate)
        t_id_dict[transcript_id]["zero_and_singleton_expectation_sum"] += float(zeroton_plus_singleton_mutation_rate)

        t_id_dict[transcript_id]["{}_zeroton_expected".format(delta_score_bin)] += float(zeroton_mutation_rate) 
        t_id_dict[transcript_id]["{}_{}_zeroton_expected".format(delta_score_bin,strand_corrected_ref)] += float(zeroton_mutation_rate) 
        t_id_dict[transcript_id]["{}_zero_and_singleton_expected".format(delta_score_bin)] += float(zeroton_plus_singleton_mutation_rate) 
        t_id_dict[transcript_id]["{}_{}_zero_and_singleton_expected".format(delta_score_bin,strand_corrected_ref)] += float(zeroton_plus_singleton_mutation_rate) 

        ## Get the best delta score bin for the current position
        per_pos_max_delta_bin[pos] = get_best_delta_bin(per_pos_max_delta_bin[pos], delta_score_bin) 

        ## Add to the positions considered 
        considered_position.add(pos)

        ## set the position as a variant added
        added_var_pos.add(pos)

    ## If no var, or var doesn't have PASS filter, or coverage is below cc, or var is non-SNV, REF->REF
    for region_pos in spliceai_region_dict.keys():
        
        ## If the region position does not overlap the query region, skip it
        if region_pos < region_start or region_pos > region_end:
            continue

        ## If the position had already been added to
        if region_pos in added_var_pos:
            continue

        ## get ref allele for current position
        ref_allele = list(spliceai_region_dict[region_pos].keys()).pop() if region_pos in spliceai_region_dict and len(spliceai_region_dict[region_pos].keys()) == 1 else "None"

        ## Get strand corrected alleles
        strand_corrected_ref = correct_allele_by_strand(region_strand, ref_allele)

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
        ### THis happends when multiple ref allels occur but no ref allele matches the fasta ref allele
        if len(spliceai_region_dict[region_pos]) == 0:
            error = ("\n\n!!ERROR!! Missing Spliceai Info"
                    " at Position = {}:{}\n"
                    " Skipping Site").format(region_chrom, region_pos)
            output_log(error,error_log_file)
            empty_position_info += 1
            cur_empty_position_info += 1
            continue
            

        ## Get the max sum of delta score among all ref to alt combinations for the current position 
        max_sum_delta_score = max([get_max_sum_delta_score(spliceai_region_dict[region_pos][ref_allele][alt_allele]) for alt_allele in spliceai_region_dict[region_pos][ref_allele].keys()])
        delta_score_bin = get_delta_score_bin(max_sum_delta_score, delta_score_bins)

        zerotons = 1.0
        zerotons_plus_singletons = 1.0
        singletons = 0.0
        non_zerotons = 0.0
        non_zerotons_plus_singletons = 0.0 

        ## Get zeroton and zeroton plus singleton mutation rate for the current delta score bin and ref allele
        zeroton_mutation_rate = float(mut_rate_dict[delta_score_bin]["zeroton"][strand_corrected_ref])
        zeroton_plus_singleton_mutation_rate = float(mut_rate_dict[delta_score_bin]["zeroton_plus_singleton"][strand_corrected_ref]) 

        ## add ovserved and expected scores
        ### Observed scores
        t_id_dict[transcript_id]["zeroton_observed_var_count"] += non_zerotons
        t_id_dict[transcript_id]["zero_and_singleton_observed_var_count"] += non_zerotons_plus_singletons

        t_id_dict[transcript_id]["{}_zeroton_observed".format(delta_score_bin)] += non_zerotons
        t_id_dict[transcript_id]["{}_{}_zeroton_observed".format(delta_score_bin,strand_corrected_ref)] += non_zerotons
        t_id_dict[transcript_id]["{}_zero_and_singleton_observed".format(delta_score_bin)] += non_zerotons_plus_singletons 
        t_id_dict[transcript_id]["{}_{}_zero_and_singleton_observed".format(delta_score_bin,strand_corrected_ref)] += non_zerotons_plus_singletons 

        ## Expected Scores
        t_id_dict[transcript_id]["zeroton_expectation_sum"] += float(zeroton_mutation_rate)
        t_id_dict[transcript_id]["zero_and_singleton_expectation_sum"] += float(zeroton_plus_singleton_mutation_rate)

        t_id_dict[transcript_id]["{}_zeroton_expected".format(delta_score_bin)] += float(zeroton_mutation_rate) 
        t_id_dict[transcript_id]["{}_{}_zeroton_expected".format(delta_score_bin, strand_corrected_ref)] += float(zeroton_mutation_rate) 
        t_id_dict[transcript_id]["{}_zero_and_singleton_expected".format(delta_score_bin)] += float(zeroton_plus_singleton_mutation_rate) 
        t_id_dict[transcript_id]["{}_{}_zero_and_singleton_expected".format(delta_score_bin,strand_corrected_ref)] += float(zeroton_plus_singleton_mutation_rate) 

        ## Get the best delta score bin for the current position
        per_pos_max_delta_bin[region_pos] = get_best_delta_bin(per_pos_max_delta_bin[region_pos], delta_score_bin) 
    
        ## Add to the positions considered 
        considered_position.add(region_pos)

    ## Update the transcript dict with error info
    t_id_dict[transcript_id]["non_matching_ref_alelles"] += cur_non_matching_ref_alleles
    t_id_dict[transcript_id]["non_matching_symbol"] += cur_non_matching_symbols
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



#---------------------------------------------------------------------------------------------------------------------------------
## Main
#---------------------------------------------------------------------------------------------------------------------------------

def gene_o_e(parser, args):


    print("\n\n\t***********************")
    print("\t* Gene O and E Scores *")
    print("\t***********************\n\n")


    start_time = time.time()
    args = arguments()

    if args.chrom == None:
        args.chrom = ["All"]

    print(("\nInput Arguments:"
           "\n================"
           "\n - gnomad-vcf:           {}"
           "\n - coverage:             {}"
           "\n - gtf-file:             {}"
           "\n - spliceai-vcf:         {}"
           "\n - alt-gene-symbol:      {}"
           "\n - fasta:                {}"
           "\n - seg-dups              {}"
           "\n - self-chains           {}" 
           "\n - out-file:             {}"
           "\n - mutation-rate-table:  {}"
           "\n - chrom:                {}"
           "\n - cc:                   {}"
           "\n - coverage-label:       {}"
           "\n - n-cpu:                {}"
           "\n - repeat-score-cutoff   {}"
           ).format(args.gnomad_vcf, 
                    args.coverage, 
                    args.gtf_file, 
                    args.spliceai_vcf, 
                    args.alt_gene_symbol,
                    args.fasta, 
                    args.seg_dups,
                    args.self_chains,
                    args.out_file, 
                    args.mutation_rate_table,
                    ", ".join(args.chrom), 
                    args.cc, 
                    args.coverage_label,
                    args.n_cpu,
                    args.repeat_score_cutoff
                    )
    )

    ## Check coverage cutoff value
    if args.cc > 1.0 or args.cc < 0.0: 
        print("\n!!ERROR!! The Coverage cutoff is set at {}. It needs to be between 0.0 and 1.0. Please correct the problem and try again".format(args.cc))
        sys.exit(1)


    ## Check number of CPUs
    ncpus = args.n_cpu if args.n_cpu > 0 else 1
    if ncpus > multiprocessing.cpu_count(): 
        print("\n**Warning** Too many CPUs designated. Number of CPUs will be set to {}".format( multiprocessing.cpu_count() if multiprocessing.cpu_count() > 0 else 1)) 
        ncpus = multiprocessing.cpu_count() if multiprocessing.cpu_count()  > 0 else 1

    
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
    gene_strand_dict = expand_gene_dict_with_alt_symbols(alt_symbol_dict, gene_strand_dict)


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
    print("\n\nGetting the mutation rates for each SpliceAI delta bin and ref allele combo")
    print("===========================================================================")
    header_index = 0
    try:
        with io.open(args.mutation_rate_table, "rt", encoding = "utf-8") as mrt:
            for i,line in enumerate(mrt):
                if line[0] == "#": 
                    header_index = i
                elif line[0] != "#":
                    break
    except IOError as e:
        print("\n!!ERROR!! Unable to read the mutation rate table")
        print(str(e))
        sys.exit(1)

    print("\n\tUsing the SpliceAI delta bin mutation count table: '{}'".format(args.mutation_rate_table))
    ## Load the nss table into a pandas df
    mr_table = pd.read_csv(args.mutation_rate_table, sep="\t", index_col=False, header = header_index)
    mr_table = mr_table.rename(columns = {"#ref":"ref"})


    print("\tGetting mutation rates")
    ## Get the mutate rate for each near splice position per ref allele
    by_ref_mutation_rate_dict = defaultdict(lambda: defaultdict(float))
    mr_table.groupby(["delta_score_bin"]).apply(get_normalized_matrix_per_label, mutation_rate_dict = by_ref_mutation_rate_dict)

    print("\nMutation rates by SpliceAI bin and reference allele")
    print("-------------------------------------------------")
    for label,value in by_ref_mutation_rate_dict.items():
        print("SpliceAI delta bin: {}".format(label))
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

    print("\n\tPreapting transcript variant tracking")
    pbar = tqdm(total = len(transcript_region_dict))
    pbar.set_description("Preparing tracking structure")
    delta_score_bin_list = list(by_ref_mutation_rate_dict.keys())
    for key in transcript_region_dict:

        pbar.update(1)

        for label in delta_score_bin_list:

            ## Add label specific observed and expected score counts
            ### set observed score to 1 as default 
            ### Get the average mutation rate for the current bin and set as default for expectation   
            transcript_region_dict[key]["{}_zeroton_observed".format(label)] = 1.0
            transcript_region_dict[key]["{}_zero_and_singleton_observed".format(label)] = 1.0
            transcript_region_dict[key]["{}_zeroton_expected".format(label)] = 1.0
            transcript_region_dict[key]["{}_zero_and_singleton_expected".format(label)] = 1.0

            ## Delta score bin count
            transcript_region_dict[key]["{}_count".format(label)] = 0.0

            ## For each label, get a count per ref allele
            for ref_allele in ["A","C","G","T"]:
                
                ## Set observed to 1
                ## Set expected to the mutation rate for that ref allele
                transcript_region_dict[key]["{}_{}_zeroton_observed".format(label,ref_allele)] = 1.0
                transcript_region_dict[key]["{}_{}_zero_and_singleton_observed".format(label,ref_allele)] = 1.0
                transcript_region_dict[key]["{}_{}_zeroton_expected".format(label,ref_allele)] = 1.0
                transcript_region_dict[key]["{}_{}_zero_and_singleton_expected".format(label,ref_allele)] = 1.0  

        ## Add general observed and expected score counts
        transcript_region_dict[key]["zeroton_observed_var_count"] = 0.0
        transcript_region_dict[key]["zero_and_singleton_observed_var_count"] = 0.0
        transcript_region_dict[key]["zeroton_expectation_sum"] = 0.0
        transcript_region_dict[key]["zero_and_singleton_expectation_sum"] = 0.0

        ## Add error tracking per transcript
        transcript_region_dict[key]["non_matching_ref_alelles"] = 0.0
        transcript_region_dict[key]["non_matching_symbol"] = 0.0
        transcript_region_dict[key]["multiple_ref_allele_pos"] = 0.0
        transcript_region_dict[key]["ref_alleles_as_N"] = 0.0
        transcript_region_dict[key]["empty_position"] = 0.0
        transcript_region_dict[key]["seg_dup_pos"] = 0.0
        transcript_region_dict[key]["self_chain_pos"] = 0.0
        transcript_region_dict[key]["positions_considered"] = 0.0
        transcript_region_dict[key]["total_positions"] = transcript_region_dict[key]["end"] - transcript_region_dict[key]["start"] + 1

    pbar.close()

    print("\n\nProcessing transcript regions") 
    print("=============================\n") 
    pbar = tqdm(total = len(transcript_region_dict))
    pbar.set_description("Processing transcript regions")
    prev_chrom = 1

    ## Track number of non matching sites
    non_matching_ref_alleles = 0
    non_matching_symbols = 0
    pos_with_multiple_refs = 0
    ref_allele_n = 0
    empty_position_info = 0
    seg_dup_or_self_chain = 0


    ## Create log files
    date_string = datetime.datetime.now().strftime("%m-%d-%Y")
    time_string = datetime.datetime.now().strftime("%H:%M:%S")
    output_log("##Start Date: {}\n##Start Time: {}\n#ERRORS:\n".format(date_string,time_string), 
               error_log_file,
               new_file = True)
              
    output_log("##Start Date: {}\n##Start Time: {}\n#STDOUT:\n".format(date_string, time_string), 
               out_log_file,
               new_file = True)

    for transcript_id in transcript_region_dict.keys():
        
        pbar.update(1)

        chrom = transcript_region_dict[transcript_id]["chrom"]
        start = transcript_region_dict[transcript_id]["start"]
        end = transcript_region_dict[transcript_id]["end"]
        strand = transcript_region_dict[transcript_id]["strand"]
        gene_name = transcript_region_dict[transcript_id]["gene_name"] 

        ## Skip any non selected chromosomes 
        if args.chrom != ["All"] and chrom not in args.chrom:
            continue

        ## reload vcf
        if prev_chrom != chrom:
            spliceai_vcf = VCF(args.spliceai_vcf, threads=ncpus)
            gnomad_vcf = VCF(args.gnomad_vcf, threads=ncpus)
            prev_chrom = chrom

        ## Get SpliceAI region Info
        spliceai_position_dict = dict()
        min_region_start = float("inf")
        max_region_start = 0
        prev_pos = -1
        spliceai_region_list = []
        multi_ref_pos = set()
        for i,spliceai_var in enumerate(spliceai_vcf("{}:{}-{}".format(chrom,start,end))):

            ## Get the spliceai annotation as a dict
            spliceai_annotation = dict(zip(spliceai_info_index, spliceai_var.INFO["SpliceAI"].strip().split("|")))

           ## Check for the correct strand
            if spliceai_annotation["SYMBOL"] in gene_strand_dict:
                if gene_strand_dict[spliceai_annotation["SYMBOL"]]["strand"] != strand:
                    ## Skip this varriant if it is on the wrong strand
                    continue
            else:
                ## If the gene name is not in the gene strand dict, where alternative symobls have also been added, skip this site
                continue

            matching_symbols = set(gene_strand_dict[gene_name]) 
            matching_symbols.add(gene_name)

            ## Skip any spliceai var that does not having a matching symobl with the current gene
            if spliceai_annotation["SYMBOL"] not in matching_symbols:
                continue
            
            ## Add chrom, pos, ref, and alt to the dict
            spliceai_annotation.update({"chrom":spliceai_var.CHROM, 
                                         "pos":spliceai_var.POS, 
                                         "ref":spliceai_var.REF, 
                                         "alt":spliceai_var.ALT[0]})

            ## Check for seg dups or self chains at the current posisiont, If overlap, skip this position
            ## If seg dup and self chain overlap, the seg dup count will be incremetned over the self chain count
            if list(segdup_interlap_dict[chrom].find((spliceai_var.POS, spliceai_var.POS))): 
                seg_dup_or_self_chain += 1
                transcript_region_dict[transcript_id]["seg_dup_pos"] +=  1.0  
                continue

            elif list(selfchain_interlap_dict[chrom].find((spliceai_var.POS, spliceai_var.POS))):
                seg_dup_or_self_chain += 1
                transcript_region_dict[transcript_id]["self_chain_pos"] += 1.0 
                continue


            ## keep track of position in the region that have SpliceAI annotations
            ## Track continous regions annotated by spliceai
            if prev_pos == -1:  
                prev_pos = spliceai_var.POS


            ## Continous region tracking
            if prev_pos == spliceai_var.POS:
                prev_pos = spliceai_var.POS
            elif prev_pos + 1 == spliceai_var.POS:
                prev_pos = spliceai_var.POS
            else:
                spliceai_region_list.append({"chrom":chrom, "start":min_region_start, "end":max_region_start})
                min_region_start = spliceai_var.POS
                max_region_start = spliceai_var.POS
                prev_pos = spliceai_var.POS

            ## relative min and max
            if spliceai_var.POS < min_region_start:
                min_region_start = spliceai_var.POS
            if spliceai_var.POS > max_region_start:
                max_region_start = spliceai_var.POS



            """
            splice position dict:
            =====================
            { Position :
                { Ref Allele :
                    { ALT Allele : 
                        { Gene Symbol :
                            { Spliceai Annotation dict }
                        }
                    }
                }
            }
            """
            ## Update spliceai position dict
            if spliceai_var.POS not in spliceai_position_dict:
                spliceai_position_dict[spliceai_var.POS] = {spliceai_var.REF:{spliceai_var.ALT[0]:{spliceai_annotation["SYMBOL"]:spliceai_annotation}}}

            else:
                ## Check for additional reference alleles. 
                if spliceai_var.REF not in spliceai_position_dict[spliceai_var.POS]:

                    ## get ref allele from fasta
                    fa = Fasta(args.fasta)
                    fa_ref_allele = fa[spliceai_var.CHROM][spliceai_var.start:spliceai_var.end]
                    fa.close()

                    ## iterate over each ref alelle for the current position
                    for ref_allele in list(spliceai_position_dict[spliceai_var.POS].keys()) + [spliceai_var.REF]:
                        ## if the referenec allele does not match the fasta reference and the key is in the dict, remove it
                        if ref_allele != fa_ref_allele and ref_allele in spliceai_position_dict[spliceai_var.POS]:
                            del spliceai_position_dict[spliceai_var.POS][ref_allele]
                        ## if the reference allele does match, but it is dont in the dict, add it
                        elif ref_allele == fa_ref_allele and ref_allele not in spliceai_position_dict[spliceai_var.POS]:
                            ## check if ref is for the current var
                            if ref_allele == spliceai_var.REF:
                                spliceai_position_dict[spliceai_var.POS] = {spliceai_var.REF:{spliceai_var.ALT[0]:{spliceai_annotation["SYMBOL"]:spliceai_annotation}}}
                        ## if ref allele matches fasta ref allele, and it is already in the dict, add alt allele
                        elif ref_allele == fa_ref_allele and ref_allele in spliceai_position_dict[spliceai_var.POS]:
                            ## check if ref is for the current var
                            if ref_allele == spliceai_var.REF:
                                if spliceai_var.ALT[0] not in spliceai_position_dict[spliceai_var.POS][spliceai_var.REF]:
                                    spliceai_position_dict[spliceai_var.POS][spliceai_var.REF].update({spliceai_var.ALT[0]:{spliceai_annotation["SYMBOL"]:spliceai_annotation}})
                                else:
                                    spliceai_position_dict[spliceai_var.POS][spliceai_var.REF][spliceai_var.ALT[0]].update({spliceai_annotation["SYMBOL"]:spliceai_annotation})
                                
                    ## Add Warning
                    error = ("\n**WARNING** SpliceAI has multiple reference alleles for a single genomic position. Checking for correct reference allele."
                            "Genomic Position: {}:{}\n").format(spliceai_var.CHROM, spliceai_var.POS)
                    if len(spliceai_position_dict[spliceai_var.POS].keys()) < 1:
                        error += ("No matching ref allele for current position."
                                  "\n\tSpliceAI Ref Alleles = {}"
                                  "\n\tFasta Ref Allele = {}").format(",".join(list(spliceai_position_dict[spliceai_var.POS].keys()) + [spliceai_var.REF]), 
                                                                      fa_ref_allele)
                    output_log(error,error_log_file)
                    if spliceai_var.POS not in multi_ref_pos:
                        pos_with_multiple_refs += 1
                        multi_ref_pos.add(spliceai_var.POS)

                ## Add alt alleles 
                elif spliceai_var.ALT[0] not in spliceai_position_dict[spliceai_var.POS][spliceai_var.REF]: 
                    spliceai_position_dict[spliceai_var.POS][spliceai_var.REF].update({spliceai_var.ALT[0]:{spliceai_annotation["SYMBOL"]:spliceai_annotation}})

                else: ## If multiple genes, add each gene separately
                    spliceai_position_dict[spliceai_var.POS][spliceai_var.REF][spliceai_var.ALT[0]].update({spliceai_annotation["SYMBOL"]:spliceai_annotation})

        ## Add the last region if it has not been added yet
        if len(spliceai_region_list) == 0 or {"chrom":chrom, "start":min_region_start, "end":max_region_start} != spliceai_region_list[-1]:
            spliceai_region_list.append({"chrom":chrom, "start":min_region_start, "end":max_region_start})

        
        ## out log
        stdout_string = ("\nBefore Filtering Region: {}:{}-{}"
                         "\nSpliceAI Filtered Region(s): {}").format(chrom,start,end,spliceai_region_list)
        output_log(stdout_string,out_log_file)

        ## the number of position with multiple ref alleles 
        transcript_region_dict[transcript_id]["multiple_ref_allele_pos"] = len(multi_ref_pos)

        ## Iterate over each relateive spliceai continous region
        for relative_spliceai_region in spliceai_region_list:

            ## Skip regions where no spliceai info was found
            if relative_spliceai_region["start"] == float("inf") and relative_spliceai_region["end"] == 0:
                continue

            ## Get a dict of coverage where key = position, value = coverage value
            ### {POS:COVERAGE,POS:COVERAGE,POS:COVERAGE,...}
            by_position_coverage_dict = {int(x.strip().split("\t")[cov_header_dict["end"]]):float(x.strip().split("\t")[cov_header_dict[args.coverage_label]])
                                            for x in cov_file.fetch("chr{}".format(relative_spliceai_region["chrom"]),
                                                                                                        int(relative_spliceai_region["start"] - 1),
                                                                                                        int(relative_spliceai_region["end"] + 1))}

            ## Update mutation frequencies from gnomAD variants
            (non_matching_ref_alleles, 
             non_matching_symbols, 
             ref_allele_n, 
             empty_position_info) = get_gnomad_vars_by_position(gnomad_vcf,
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
                                                                transcript_region_dict,
                                                                transcript_id,
                                                                by_ref_mutation_rate_dict,
                                                                non_matching_ref_alleles,
                                                                non_matching_symbols,
                                                                ref_allele_n,
                                                                empty_position_info)
        

    ## Close progress bar
    pbar.close()

    end_time = time.time()

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




    ## Write mutation table
    print("\nWriting by transcript observed and expected scores to '{}'".format(args.out_file))

    with open(args.out_file, "w") as out: 
        out_header = ["#chrom",
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

    
        delta_score_bin_list = list(by_ref_mutation_rate_dict.keys())
        out_header.extend(["{}_zeroton_observed".format(label) for label in delta_score_bin_list])
        out_header.extend(["{}_zeroton_expected".format(label) for label in delta_score_bin_list])
        out_header.extend(["{}_zero_and_singleton_observed".format(label) for label in delta_score_bin_list])
        out_header.extend(["{}_zero_and_singleton_expected".format(label) for label in delta_score_bin_list])
        out_header.extend(["{}_count".format(label) for label in delta_score_bin_list])

        for d_bin in delta_score_bin_list:
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
        out_header.append("seg_dup_pos")
        out_header.append("self_chain_pos")
        out_header.append("positions_considered")
        out_header.append("total_positions")
        out.write("\t".join(out_header) + "\n")

        
        for key,value in transcript_region_dict.items():

            output_list = [str(value["chrom"]), 
                           str(value["start"]), 
                           str(value["end"]), 
                           str(value["strand"]), 
                           str(value["gene_id"]), 
                           str(value["gene_name"]), 
                           str(key),
                           str(value["max_exon_number"]),
                           str(value["cds_exon_count"]),
                           str(value["zeroton_observed_var_count"]), 
                           str(value["zeroton_expectation_sum"]),
                           str(value["zero_and_singleton_observed_var_count"]),
                           str(value["zero_and_singleton_expectation_sum"])
                           ]    

            output_list.extend([str(value["{}_zeroton_observed".format(label)]) for label in delta_score_bin_list])
            output_list.extend([str(value["{}_zeroton_expected".format(label)]) for label in delta_score_bin_list])
            output_list.extend([str(value["{}_zero_and_singleton_observed".format(label)]) for label in delta_score_bin_list])
            output_list.extend([str(value["{}_zero_and_singleton_expected".format(label)]) for label in delta_score_bin_list])
            output_list.extend([str(value["{}_count".format(label)]) for label in delta_score_bin_list])

            for d_bin in delta_score_bin_list:
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
            output_list.append(str(value["seg_dup_pos"]))
            output_list.append(str(value["self_chain_pos"]))
            output_list.append(str(value["positions_considered"]))
            output_list.append(str(value["total_positions"]))

            out.write("\t".join(output_list) + "\n")
            output_list = []

    print("\nDONE")

