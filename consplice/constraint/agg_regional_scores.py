import sys
import os
import pandas as pd
import numpy as np
import argparse 
from scipy.stats.mstats import gmean
from scipy.stats import hmean
from tqdm import tqdm
from interlap import InterLap
from collections import defaultdict
from .calculate_constraint_score import convert_scores_to_percentiles


#---------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
#---------------------------------------------------------------------------------------------------------------------------------

def add_aggregate_overlapping_regions(sub_p):
    

    p = sub_p.add_parser("agg-overlapping-reg",
                         help = "Aggregate scores from overlapping regions where the step size of a region is less than the window size of the region",
                         formatter_class=argparse.RawDescriptionHelpFormatter,
                         description=("\n\t***************************************\n"
                                      "\t* Aggregate Overlapping Region Scores *\n"
                                      "\t***************************************\n\n"
                                      "\tCombine overlapping regional constraint scores into an aggergated score\n"
                                      "\tusing the region's step size. All regions that intersect the user defined\n"
                                      "\tstep size will be combined together into a single score."
                        )

    )

    req = p.add_argument_group("Required Arguments")
    req.add_argument(
        "--score-file",
        metavar = "Regional Score File", 
        required=True, 
        help="(Required) The path to the regional score file to combine scores for"
    )

    req.add_argument(
        "--score-field",
        metavar = "O/E Score Field",
        required = True,
        help = "(Required) The name of the O/E score field within the regional score file to combine based on intersecting regions."
    )

    req.add_argument(
        "--step-size",
        metavar = "Step Size",
        required=True,
        help="(Required) The step size used to create the overlapping regions. This step size is used to step through each gene in the score file and the scores from all region that intersects with a step region will be combined together."
    )   

    req.add_argument(
        "--out-file",
        metavar = "Output file",
        required = True,
        help = "(Required) The name of the output file to create."
    )

    p.add_argument(
        "--new-score-name",
        metavar="New O/E Score Column Name",
        default = "Aggregated_ConSplice_O/E",
        help="(Optional) The name of the new score column that is created. Default = 'Aggregated_ConSplice_O/E'"
    )

    p.add_argument(
        "--new-pctl-name",
        metavar="New Percentile Score Column Name",
        default = "Aggergated_ConSplice_Percentile",
        help="(Optional) The name of the new percentile score column that is created. Default = 'Aggergated_ConSplice_Percentile'"
    )

    p.add_argument(
        "--invert-pctl",
        action="store_true",
        help = "(Optional) Whether or not to invert the new scores before converting the scores to percentiles. By default, the percentile score will correlate with increasing scores from the new aggregated score created during the intersection process. If the smaller (and/or more negative) a score is correlates with a higher percentile then this argument should be set"

    )

    p.add_argument(
        "--agg-type",
        choices = ["mean","median"],
        default = "median",
        help = "(Optional) How to aggregate the overlapping scores together. Choices = 'mean' or 'median'. Default = 'median'"
    )


    p.set_defaults(func=agg_overlapping_regions)


#---------------------------------------------------------------------------------------------------------------------------------
## Functions/Methods
#---------------------------------------------------------------------------------------------------------------------------------

def create_interlap_object(df, 
                           chrom_col, 
                           start_col, 
                           end_col, 
                           tx_start_col, 
                           tx_end_col, 
                           strand_col,
                           tx_id_col,
                           gene_id_col,
                           gene_symbol_col,
                           score_col):
    """
    create_interlap_object
    ======================
    Method to create an interlap object, interval tree with little overhead. This interlap object will take start and end positions 
     from a pandas dataframe and create an interlap object for each interval. Interlap objects will be seperated into chromosomes using 
     the chromosome designation in the dataframe. Additional gene information for each region will be added. For more information about
     InterLap and how to use it see: https://brentp.github.io/interlap/

    Parameters:
    -----------
    1) df: (Pandas DataFrame) A pandas dataframe where each row represents an interval. This DF should have the following columns: chrom, start, end, transcript start, transcript end, 
                               strand, transcript id, gene id, gene symbol, and score column
    2) chrom_col:       (str) The name of the column in the dataframe that has the chromosome name
    3) start_col:       (str) The name of the column in the dataframe that has represents the start positions of the interval
    3) end_col:         (str) The name of the column in the dataframe that has represents the end positions of the interval
    4) tx_start_col:    (str) The name of the column in the dataframe that has represents the start positions of the transcipt
    5) tx_end_col:      (str) The name of the column in the dataframe that has represents the end positions of the transcipt
    6) strand_col:      (str) The name of the column in the dataframe that has represents the strand. (+ or -)
    7) tx_id_col:       (str) The name of the column in the dataframe that has represents the transcript ID 
    8) gene_id_col:     (str) The name of the column in the dataframe that has represents the gene ID 
    9) gene_symbol_col: (str) The name of the column in the dataframe that has represents the gene symbol 
    10) score_col:      (str) The name of the column in the dataframe that has represents the desired score to keep  

    Returns:
    ++++++++
    1) (dict) A dictionary of interlap objects. Keys = chromsome names. Values = Interlap object created from intervals in the dataframe with the associated chromosome. 
    """
    
    ## Create a dictionary that will constraint interlap objects per chromosome
    interlap_dict = defaultdict(InterLap)
    
    ## Iterate over the dataframe by row
    for row in df.itertuples():

        ## Add row contents to the chromosome specific interlap object
        interlap_dict[getattr(row,chrom_col)].add((int(getattr(row,start_col)),
                                                   int(getattr(row,end_col)), 
                                                   getattr(row,tx_start_col), 
                                                   getattr(row,tx_end_col), 
                                                   getattr(row,strand_col), 
                                                   getattr(row,tx_id_col), 
                                                   getattr(row,gene_id_col), 
                                                   getattr(row,gene_symbol_col), 
                                                   getattr(row,score_col)))

    return(interlap_dict)
        

def aggregate_intersecting_regions(chrom, start_pos, end_pos, strand, symbol, step_size, interlap_dict, agg_type):
    """
    aggregate_intersecting_regions
    ==============================
    Method to get the overlapping/intersecting interval from an interlap object and different regions in a gene and aggegate the scores for 
     all overlapping interval. A new region will be created within a gene using a step size, the region will be checked for intersecting intervals 
     in the interlap object, and the scores from all overlapping intervals will be combined into one score for that region. This is done across the 
     gene/region using the step size until the entire gene/region has scores for each step.

    Parameters:
    -----------
    1) chrom:         (str)  The name of the chromosome of the gene/region of interest
    2) start_pos:     (int)  The start position of the region or gene
    3) end_pos:       (int)  The end position of the region or gene
    4) strand:        (str)  The strand of the region or gene (+ or -)
    5) symbol:        (str)  The gene symbol for the region or gene
    6) step_size:     (int)  The size of regions to create accross the gene/region
    7) interlap_dict: (dict) A dictionary of interlap objects created from the 'create_interlap_object' method
    8) agg_type       (str)  The aggregation type. (mean or median)

    Returns:
    ++++++++
    1) (list) A 2d list, where each inner list represents a new region with associated information in the following order:
                chromosome, start, end, transcript start, transcript end, strand, transcript id, gene id, gene symbol, region size, aggregated score
    """

    new_regions = []

    temp_start = start_pos
    temp_end = start_pos + (step_size - 1)

    ## Get each intersecting regions for each step 
    while temp_end < end_pos:

        o_over_e_values = []
        tx_start = 0
        tx_end = 0
        tx_id = ""
        gene_id = ""

        for intersect in interlap_dict[chrom].find((temp_start,temp_end)):
            ## Check that the intersect object:
            ## 1) the strand of the intersect object is the same as the region of interset 
            ## 2) the gene symbols match
            ## 3) Does not have a matching end position to the start position of the step. (These should not be included)
            if intersect[4] == strand and intersect[7] == symbol and intersect[1] != temp_start:
                o_over_e_values.append(intersect[8])
                tx_start = intersect[2] if tx_start == 0 else tx_start
                tx_end = intersect[3] if tx_end == 0 else tx_end
                tx_id = intersect[5] if tx_id == "" else tx_id
                gene_id = intersect[6] if gene_id == "" else gene_id

        if o_over_e_values:
            new_regions.append([chrom, temp_start, temp_end, tx_start, tx_end, strand, tx_id, gene_id, symbol, step_size, np.mean(o_over_e_values) if agg_type == "mean" else np.median(o_over_e_values)]) 

        temp_start += step_size
        temp_end = temp_start + (step_size - 1)

    ## Get the last region based on the end pos
    ## Genes that are smaller then the step size will not have any new regions 
    ### Example of some genes with small gene sizes (~less than 300bp) SPRR1A, LCE3A, LCE3B, LCE1F, MT1HL1, PYDC2, KRTAP22-1, KRTAP19-6, etc.
    if new_regions and new_regions[-1][2] < end_pos:
        temp_start = new_regions[-1][2] + 1
        temp_end = end_pos
        new_end = 0
        for intersect in interlap_dict[chrom].find((temp_start,temp_end)):
            if intersect[4] == strand and intersect[7] == symbol:

                tx_start = intersect[2] if tx_start == 0 else tx_start
                tx_end = intersect[3] if tx_end == 0 else tx_end
                tx_id = intersect[5] if tx_id == "" else tx_id
                gene_id = intersect[6] if gene_id == "" else gene_id

                ## Set the end of the new region to the end of the longest overlapping region
                if intersect[1] > new_end:
                    new_end = intersect[1]

        if o_over_e_values:
            new_regions.append([chrom, temp_start, new_end, tx_start, tx_end, strand, tx_id, gene_id, symbol, step_size, np.mean(o_over_e_values) if agg_type == "mean" else np.median(o_over_e_values)]) 

    return(new_regions)


#---------------------------------------------------------------------------------------------------------------------------------
## Main
#---------------------------------------------------------------------------------------------------------------------------------

def agg_overlapping_regions(parser, args):


    print("\n\t***************************************")
    print("\t* Aggregate Overlapping Region Scores *")
    print("\t***************************************\n\n")


    print(("\nInput Arguments:"
           "\n================"
           "\n - score-file:           {}"
           "\n - score-field:          {}"
           "\n - step-size:            {}"
           "\n - out-file:             {}"
           "\n - new-score-name:       {}"
           "\n - invert-pctl:          {}" 
           "\n - new-pctl-score-name:  {}"
           "\n - agg-type:             {}"
           ).format(args.score_file, 
                    args.score_field, 
                    args.step_size, 
                    args.out_file, 
                    args.new_score_name,
                    args.invert_pctl,
                    args.new_pctl_name,
                    args.agg_type,
                    )
    )

    print("\n> Reading score file")
    ## Read scores file into dataframe
    scores_df = pd.read_csv(args.score_file, sep = "\t", index_col = False)

    ## Change chrom column
    scores_df =  scores_df.rename(columns = {"#chrom":"Chromosome"})

    print("\n> Creating InterLap Object")
    ## Create a dictionary where keys are chromosomes and the value is an interlap object for regions in that chromosome
    ### These interlap objects will be used for quick interval intersection look up
    inter_dict = create_interlap_object(scores_df, 
                                        "Chromosome",
                                        "region_start",
                                        "region_end",
                                        "txStart",
                                        "txEnd",
                                        "strand",
                                        "transcript_id",
                                        "gene_id",
                                        "gene_symbol",
                                        args.score_field,
                                        )

    ## get the number of genes
    ngenes = scores_df.groupby(["gene_id","transcript_id"]).first().shape[0]

    print("\n> Aggregating regions by step size: {}".format(args.step_size))

    ## Iterate over each gene and get the new step size regions with the aggergated score
    step_regions = []
    for i,row in enumerate(scores_df.groupby(["gene_id","transcript_id"]).first().itertuples()):

        step_regions.extend(aggregate_intersecting_regions(chrom = row.Chromosome, start_pos = int(row.txStart), end_pos = int(row.txEnd), strand = row.strand, symbol = row.gene_symbol, step_size = int(args.step_size), interlap_dict = inter_dict), agg_type = args.agg_type)

    ## Create new dataframe
    print("\n> Combining all regions")
    new_regions_df = pd.DataFrame(step_regions, columns = ["#chrom","region_start","region_end","txStart","txEnd","strand","transcript_id","gene_id","gene_symbol","step_size",args.new_score_name])

    print("\n> Converting new scores to percentiles")
    new_region_w_pctl_df = convert_scores_to_percentiles(query_df = new_regions_df, 
                                                  score_column = args.new_score_name, 
                                                  percentile_column_name = args.new_pctl_name, 
                                                  invert_percentiles = args.invert_pctl) 

    print("\n> Sorting region by alphanumeric chromosome and start position")
    new_region_w_pctl_df = new_region_w_pctl_df.sort_values(by = ["#chrom", "region_start"])

    print("\n> Wriring output to '{}'".format(args.out_file))
    new_region_w_pctl_df.to_csv(args.out_file, sep="\t", index=False)

    print("\n> DONE")

