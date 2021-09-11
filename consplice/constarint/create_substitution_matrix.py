from __future__ import print_function
import sys
import io
import argparse 
import gzip
import copy
from cyvcf2 import VCF 
import pysam
import pandas as pd
import numpy as np
import pyranges as pr
from collections import defaultdict
from tqdm import tqdm
import multiprocessing
import datetime
from pyfaidx import Fasta
import time
from .utils import (create_interlap_from_ggd_pkg,
                   correct_allele_by_strand, 
                   get_alternative_gene_symbols, 
                   get_delta_score_bin, 
                   get_max_delta_score, 
                   get_max_sum_delta_score,
                   output_log)


#---------------------------------------------------------------------------------------------------------------------------------
## Global Vars
#---------------------------------------------------------------------------------------------------------------------------------
## SpliceAI delta bins
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
error_log_file = file_date + "_count_matrix.ERROR.log"
out_log_file = file_date + "_count_matrix.OUT.log"

#---------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
#---------------------------------------------------------------------------------------------------------------------------------

def add_substitution_matrix(sub_p):

    p = sub_p.add_parser("sub-matrix",
                         help = "Build a substitution matrix using gnomAD variants and SpliceAI scores",
                         formatter_class=argparse.RawDescriptionHelpFormatter,
                         description = ("\n\t***********************\n\t* Substitution Matrix *\n\t***********************\n\n"
                                        "\tCreate a splicing-aware substitution matrix"
                                        "\n\tfrom variation in gnomAD combined with per-"
                                        "\n\tnucleotide splicing predictions from SpliceAI."
                                       ),
    )


    req = p.add_argument_group("Required Arguments")
    req.add_argument(
        "--gnomad-vcf",
        metavar = "gnomAD VCF/BCF File", 
        required=True, 
        help="(Required) The path to and name of the gnomAD vcf/bcf file."
    )

    req.add_argument(
        "--coverage",
        metavar = "gnomAD Coverage file",
        required = True,
        help = "(Required) The gnomAD coverage file for the gnomAD vcf file. The coverage files should be in a bed formated, bgzipped, and tabixed."
    )

    req.add_argument(
        "--gtf-file",
        metavar = "GTF file",
        required=True,
        help="(Required) A gtf file to get gene specific genomic coordinance from to build the count matrix."
    )   

    req.add_argument(
        "--spliceai-vcf",
        metavar = "SpliceAI SNV prediction File", 
        required=True, 
        help="(Required) The path to and name of the SpliceAI delta score predictions for every snv (vcf file)."
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
        help = "A number between 0.0 and 1.0 that represents the required number of samples with coverage at each site looked at in gnomAD. Default = 0.5."
    )

    p.add_argument(
        "--coverage-label",
        metavar = "Coverage Label",
        default = "over_10",
        help = "The label/column in the coverage file to use for the selected position coverage. Default = 'over_10'."
    )

    p.add_argument(
        "--n-cpu",
        metavar = "Number of CPUs",
        default = 3,
        type = int,
        help = "The number of CPUs to use for multi-threading. (Default = 3)."
    )

    p.add_argument(
        "--repeat-score-cutoff",
        metavar = "Repeat Score Cutoff",
        default = 0.95,
        type = float,
        help = "The cutoff score to use for seg dups or self chain repeats. Regions of a gene with a seg dup or self chain at or above this cutoff will be skipped while regions of a gene below this cutoff will not be skipped during O and E score calculations. (Default = 0.95, meaning 95 PERCENT. The score will be 0.95 for seg dups and 95.0 for self chains)."
    )

    p.set_defaults(func=substitution_matrix)


#---------------------------------------------------------------------------------------------------------------------------------
## Functions/Methods
#---------------------------------------------------------------------------------------------------------------------------------

def expand_gene_dict_with_alt_symbols(alt_symbol_dict, gene_strand_dict):
    """
    expand_gene_dict_with_alt_symbols
    =================================
    Method to update a dictionary of gene names and alternative symbols with gene strand 
    Parameters:
    -----------
    1) alt_symbol_dict:  (dict) A dictionary of possible gene symbols/names with associated symbols/names for the same gene
    2) gene_strand_dict: (dict) A dictionary of gene names with stand information
    Returns:
    ++++++++
    1) (dictionary) An updated dictionary with all alternative symbols assigned a strand if possible
    """
    
    ## Iterate over all genes in the gene strand dict
    for key in list(gene_strand_dict.keys()) :
        
        ## check if the gene exists in the alt symobl dict
        if key in alt_symbol_dict:
            
            ## Check all alt symobls 
            for gene_name in alt_symbol_dict[key]:
                
                ## If alt symobl is not in the gene strand dict, add it
                if gene_name not in gene_strand_dict:
                    
                    gene_strand_dict[gene_name] = {"gene_name":gene_name, 
                                                   "gene_id":gene_strand_dict[key]["gene_id"], 
                                                   "Strand": gene_strand_dict[key]["Strand"]}

    return(gene_strand_dict)

            
def get_merged_gtf_regions(gtf_file, feature="transcript", biotype = "protein_coding", cpus=3):
    """
    get_merged_gtf_regions
    ======================
    Method to get regions from a gtf file based on feature type and biotype of that feature. By default, 
     all regions that are protien-coding transcripts will be returned. All overlapping regions are merged into a 
     single region to remove double counting a region more than once. 

    Parameters:
    -----------
    1) gtf_file: (str) The path to and name of the gtf file to extract feature regions from 
    2) feature:  (str) The feature to extract from the gtf file (Default = transcript)
    3) biotype:  (str) The gene's biotpye to filter for. Only features with this gene biotype will be returend. (Default = protein_coding)
    4) cpu:      (int) The number of cpus to use during gtf parsing. 

    Returns:
    +++++++
    1) (pandas DataFrame) A dataframe of the gtf regions that passed filtering.    
    2) (dictionary) A dictionary of with key as the gene name of kept genes, and values as the gene name and the strand  
    """
    
    ## Load gtf into pyranges object
    gtf_pr = pr.read_gtf(gtf_file)

    ## get the biotype column
    gtf_pr_columns = gtf_pr.columns.tolist()
    biotype_column = "gene_biotype" if "gene_biotype" in gtf_pr_columns else "gene_type" if "gene_type" in gtf_pr_columns else None 

    ## Filter gtf by gene feature and biotype to get a dict of genes with strand info
    feature_gene_gtf_pr = gtf_pr.subset(lambda x: (x.Feature == "gene") & (x[biotype_column] == biotype))
    t_gene_df = feature_gene_gtf_pr.df[["gene_name","gene_id","Strand"]].transpose()
    t_gene_df.columns = t_gene_df.iloc[0]

    ## key = gene symbol / gene_name, value = dict( "gene_name": gene_name, "Strand": + or -)
    gene_name_dict = t_gene_df.to_dict() 

    ## Delete extra data 
    del feature_gene_gtf_pr
    del t_gene_df 
    
    ## Filter gtf by feature and biotype 
    filtered_gtf_pr = gtf_pr.subset(lambda x: (x.Feature == feature) & (x[biotype_column] == biotype))

    ## merge overlaping feature regions based on strand. That is, only merge regions that are on the same strand
    merged_pr = filtered_gtf_pr.merge(count = True, count_col = "Merged_Features", strand = True)

    ## Delete extra data 
    del filtered_gtf_pr 
    del gtf_pr

    ## Return the merged pyranges object as a pandas dataframe 
    return(merged_pr.df, gene_name_dict)


def create_mut_freq_table(d_score_bins):
    """
    create_mut_freq_table
    =====================
    Method to create a mutation frequency table for each of the SpliceAI delta bins and different ref/alt combinations.
     This is a raw table with default values. This table will be added to as the mutation observations occure.

    Parameters:
    -----------
    1) d_score_bins: (list) A list of 'bins' or range values to use for the new table.

    Returns:
    ++++++++
    1) (dictionary) A dictionary that represent mutation rate tracking based on different ref/alt combination and different  bins

        dict:
            { Ref Allele:
                { Alt Allele:
                    {Delta Score Bin:
                        { AC value dict}
                    }
                }
            }
    """
    
    ## Get a list of AC count dicst for the number of uniqe splice region labels
    ac_type_list = {"AC<1":0,
                    "AC<2":0,
                    "AC=1":0,
                    "AC>0":0,
                    "AC>1":0} 

    alleles = ["A","C","G","T"]
    mutation_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))

    ## Initate dict
    for ref_to_alt in [[x,y] for  x in alleles for y in alleles]:

        mutation_dict[ref_to_alt[0]][ref_to_alt[1]] = {delta_score_bin : copy.deepcopy(ac_type_list) for delta_score_bin in d_score_bins}  


    return(mutation_dict)


def get_gnomad_vars_by_position(vcf,
                                contains_chr,
                                region_chrom,
                                region_start,
                                region_end,
                                region_strand,
                                by_pos_cov_dict,
                                cov_file,
                                cov_header,
                                coverage_label,
                                coverage_cutoff,
                                vep_index,
                                spliceai_region_dict,
                                alt_symbol_dict,
                                mutation_count_table,
                                non_matching_ref_alleles,
                                non_matching_symbols,
                                ref_allele_n,
                                empty_position_info):
    """
    get_gnomad_vars_by_position
    ===========================

    Method to genomad variants by postiion and update the mutation frequency table. 
    It takes a cyvcf2 object and iterates over each varaint that is found between the position interval. Checks are done 
     for each varaint found to ensure the variant is appropriate. Checks include: 
        1) If the variant is an SNV. If the variant is not an SNV the variant will be skipped. 
        2) If the variant has a PASS filter. If the variant does not have a PASS filter the variant is skipped
        3) If the variant reference allele matches the reference allele anotated in the near splice site being compared
        4) If the variant has sufficient sequencing coverage. If not, the variant is skipped.
        5) If the variant has a matching annotated symbol. If not, varaint is skipped
    
    If a variant passes all filters it will be added to the mutation frequency table under the appropriate ref to alt keys, 
     and the splice delta score bin. 
     Example of Ref to Alt key: A->T
     Example of delta score bin: 0.8-1.0

    Info added include:
        zerotons count: (AC < 1) The number of sites without a variant 
        zerotons_plus_singletons count: (AC < 2) The number of sites with either a rare singleton variant or no variant 
        singeltons count: (AC == 1) The number of sites with a singelton variant 
        non_zerotons counts: (AC > 0) The number of sites with a variant 
        non_zerotons_plus_singletons count: (AC >1) The number of sites with a variant excluding singleton variants

    Parameters:
    -----------
    1)  vcf:                      (cyvcf2 VCF object) A cyvcf2 object for the current vcf/bcf file being queried 
    2)  contains_chr:             (bool) whether or not the the sequence names in the vcf start with chr or not
    3)  region_chrom:             (str)  The chromsome to look in for the position interval 
    4)  region_start:             (str)  The start position of the interval
    5)  region_end:               (str)  The end position of the interval
    6)  region_strand:            (str)  The strand for the current region
    7)  by_pos_cov_dict:          (dict) A dictionary of var coverage for the current region
    8)  cov_file:                 (Pysam Tabix Object) The coverage file loaded as a pysam object`
    9)  cov_header:               (list) A list representing the coverage file header
    10) coverage_label:           (str) The label that represents the coverage column in the coverage file
    11) coverage_cutoff:          (float) The coverage value at each position to use as a cutoff from the coverage file. (Between 0.0 and 1.0) 
    12) vep_index:                (list) A list that represents the vep csq index in the gnomad vcf 
    13) spliceai_region_dict:     (dict) A dictionary of SpliceAI information for the current region 
    14) alt_symbol_dict:          (dict) A dictionary of different alternative gene symbols 
    15) mutation_count_table:     (dict) The mutation frequency table for each ref to alt allele at each SpliceAI delta bin.  
    16) non_matching_ref_alleles: (int)  The current count of non matching ref alleles
    17) non_matching_symobls:     (int)  The current count of non matching symbols
    18) ref_allele_n:             (int)  The current count of reference alleles that are 'N'
    19) empty_position_info:      (int)  The current count of positions with null information 


    Returns:
    ++++++++
    1) (dir) Updated mutation_count_table,
    2) (int) Updated non_matching_ref_alleles 
    3) (int) Updated non_matching_symbols 
    4) (int) Updated ref_allele_n
    5) (int) Update empty_position_info

    """

    ## Iterate over varaints that are with the same chromosome and position
    added_var_pos = set()
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
            error = ("\n\n**WARNING** NONE matching gene symbols."
                     "\nPOSITION: {}"
                     "\n\tSite Gene: {}"
                     "\n\tVariant Gene: {}"
                     "\nSite will be skipped.").format(pos, ", ".join(spliceai_region_dict[pos][ref][alt].keys()), ", ".join(var_symbols))
            output_log(error,error_log_file)
            non_matching_symbols += 1

            added_var_pos.add(pos)
            continue

        ## Get strand corrected alleles
        strand_corrected_ref = correct_allele_by_strand(region_strand, ref)
        strand_corrected_alt = correct_allele_by_strand(region_strand, alt)

        ## Get the sum of delta scores for the current position, reference allele, and alternative allele

        max_sum_delta_score = get_max_sum_delta_score(spliceai_region_dict[pos][ref][alt])
        delta_score_bin = get_delta_score_bin(max_sum_delta_score, delta_score_bins)


        ## AC type
        zerotons = 0
        zerotons_plus_singletons = 1 if AC <= 1 else 0
        singletons = 1 if AC == 1 else 0
        non_zerotons = 1 if AC > 0 else 0
        non_zerotons_plus_singletons = 1 if AC > 1 else 0 

        ## Add count to dict
        mutation_count_table[strand_corrected_ref][strand_corrected_alt][delta_score_bin]["AC<1"] += zerotons
        mutation_count_table[strand_corrected_ref][strand_corrected_alt][delta_score_bin]["AC<2"] += zerotons_plus_singletons
        mutation_count_table[strand_corrected_ref][strand_corrected_alt][delta_score_bin]["AC=1"] += singletons
        mutation_count_table[strand_corrected_ref][strand_corrected_alt][delta_score_bin]["AC>0"] += non_zerotons
        mutation_count_table[strand_corrected_ref][strand_corrected_alt][delta_score_bin]["AC>1"] += non_zerotons_plus_singletons

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
            continue

        ## skip sites with no info
        ### THis happends when multiple ref allels occur but no ref allele matches the fasta ref allele
        if len(spliceai_region_dict[region_pos]) == 0:
            error = ("\n\n!!ERROR!! Missing Spliceai Info"
                    " at Position = {}:{}\n"
                    " Skipping Site").format(region_chrom, region_pos)
            output_log(error,error_log_file)
            empty_position_info += 1
            continue
            
        ## Get the max sum of delta scores among all ref to alt combinations for the current position 
        max_sum_delta_score = max([get_max_sum_delta_score(spliceai_region_dict[region_pos][ref_allele][alt_allele]) for alt_allele in spliceai_region_dict[region_pos][ref_allele].keys()])
        delta_score_bin = get_delta_score_bin(max_sum_delta_score, delta_score_bins)

        zerotons = 1
        zerotons_plus_singletons = 1
        singletons = 0
        non_zerotons = 0
        non_zerotons_plus_singletons = 0 

        ## Add count to dict
        mutation_count_table[strand_corrected_ref][strand_corrected_ref][delta_score_bin]["AC<1"] += zerotons
        mutation_count_table[strand_corrected_ref][strand_corrected_ref][delta_score_bin]["AC<2"] += zerotons_plus_singletons
        mutation_count_table[strand_corrected_ref][strand_corrected_ref][delta_score_bin]["AC=1"] += singletons
        mutation_count_table[strand_corrected_ref][strand_corrected_ref][delta_score_bin]["AC>0"] += non_zerotons
        mutation_count_table[strand_corrected_ref][strand_corrected_ref][delta_score_bin]["AC>1"] += non_zerotons_plus_singletons

    return(mutation_count_table, 
           non_matching_ref_alleles, 
           non_matching_symbols, 
           ref_allele_n, 
           empty_position_info)


def write_mutation_count_table(mutation_count_table, output_file):
    """
    write_mutation_count_table
    ==========================
    Method used to write the mutation frequency table to a file

    Parameters:
    -----------
    1) mutation_count_table: (dict) The count table for varaints per ref/alt/delta bin 
    2) output_file:          (str) The output file to write to
    """

    with open(output_file, "w") as out:
        ## Write header
        out.write("##zerotons = n(AC < 1)\n")
        out.write("##zerotons_plut_singletons = n(AC < 2)\n") 
        out.write("##singletons = n(AC == 1)\n")
        out.write("##non_zerotons = n(AC > 0)\n")
        out.write("##non_zerotons_plus_singletons = n(AC > 1)\n")
        out.write("#ref\talt\tdelta_score_bin\tzerotons\tzerotons_plus_singletons\tsingletons\tnon_zerotons\tnon_zerotons_plus_singletons\n")
        
        ## Write data
        for ref, alt_dict in mutation_count_table.items():
            
            for alt, delta_bin_dict in alt_dict.items():
                
                for delta_bin, ac_dict_values in delta_bin_dict.items():
                
                    out.write("\t".join([ref,
                                         alt,
                                         delta_bin,
                                         str(ac_dict_values["AC<1"]),
                                         str(ac_dict_values["AC<2"]),
                                         str(ac_dict_values["AC=1"]),
                                         str(ac_dict_values["AC>0"]),
                                         str(ac_dict_values["AC>1"])]) + "\n")


#---------------------------------------------------------------------------------------------------------------------------------
## Main
#---------------------------------------------------------------------------------------------------------------------------------

def substitution_matrix(parser, args):

    print("\n\n\t***********************************")
    print("\t* ConSplice - Substitution Matrix *")
    print("\t***********************************\n\n")

    start_time = time.time()

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

    
    print("\nGathering all strand based unique transcript regions")
    print("====================================================")
    ## Get all unique transcript regions for a strand.
    ## Additionly, get a dict of gene symbols as keys and a dict with "Strand" info as a value key
    regions_df, gene_strand_dict  = get_merged_gtf_regions(args.gtf_file,
                                                           feature="transcript",
                                                           biotype ="protein_coding", 
                                                           cpus=ncpus)
    print("\nQuery Regions:")
    print(regions_df)


    print("\nCreating an map of alternative gene symbols")
    print("===========================================")
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

    print("\nLoading vcf files")
    print("=================")

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

    print("\nPreparing processing")
    print("====================")

    print("\n\tCreating frequency table")
    mutation_frequency_dict = create_mut_freq_table(delta_score_bins) 


    print("\nProcessing each unique gene region") 
    print("==================================") 
    pbar = tqdm(total = regions_df.shape[0])
    pbar.set_description("Processing unique regions")
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

    for index,region in enumerate(regions_df.itertuples()):
        
        pbar.update(1)

        chrom = region.Chromosome
        start = region.Start
        end = region.End

        ## Skip any non selected chromosomes 
        if args.chrom != ["All"] and chrom not in args.chrom: 
            continue

        ## reload vcf to controll memmory usage
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
                if gene_strand_dict[spliceai_annotation["SYMBOL"]]["Strand"] != region.Strand:
                    ## Skip this varriant if it is on the wrong strand
                    continue
            else:
                ## If the gene name is not in the gene strand dict, where alternative symobls have also been added, skip this site
                continue

            
            ## Add chrom, pos, ref, and alt to the dict
            spliceai_annotation.update({"chrom":spliceai_var.CHROM, 
                                         "pos":spliceai_var.POS, 
                                         "ref":spliceai_var.REF, 
                                         "alt":spliceai_var.ALT[0]})


            ## Check for seg dups or self chains at the current posisiont, If overlap, skip this position
            if (list(segdup_interlap_dict[chrom].find((spliceai_var.POS, spliceai_var.POS))) 
                 or 
                 list(selfchain_interlap_dict[chrom].find((spliceai_var.POS, spliceai_var.POS)))
            ):
                seg_dup_or_self_chain += 1
                continue


            ## keep track of position in the region that have SpliceAI annotations
            ## Track continous regions annotated by spliceai
            if prev_pos == -1:  
                prev_pos = spliceai_var.POS

            ## Continous region tracking
            if prev_pos == spliceai_var.POS:
                pass
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
        stdout_string = ("\nBefore Filtering Region: {}"
                         "\nSpliceAI Filtered Region(s): {}").format(region,spliceai_region_list)
        output_log(stdout_string,out_log_file)



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
            (mutation_frequency_dict, 
             non_matching_ref_alleles, 
             non_matching_symbols, 
             ref_allele_n, 
             empty_position_info) = get_gnomad_vars_by_position(gnomad_vcf,
                                                                contains_chr, 
                                                                relative_spliceai_region["chrom"],
                                                                relative_spliceai_region["start"],
                                                                relative_spliceai_region["end"],
                                                                region.Strand,
                                                                by_position_coverage_dict,
                                                                cov_file,
                                                                cov_header,
                                                                args.coverage_label,
                                                                args.cc,
                                                                gnomad_vep_index,
                                                                spliceai_position_dict,
                                                                alt_symbol_dict,
                                                                mutation_frequency_dict,
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
    print("\nWriting mutation count table by reference allele to '{}'".format(args.out_file))
    write_mutation_count_table(mutation_frequency_dict, args.out_file)
    print("\nDONE")
        
