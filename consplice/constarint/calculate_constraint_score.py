from __future__ import print_function
import sys
import io
import os
import copy
import argparse 
import pandas as pd
import numpy as np
from scipy import stats
from collections import defaultdict



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

by_ref_delta_score_bins = ["{}_{}".format(x,y) for x in delta_score_bins for y in ["A","C","G","T"]]

recovery_choices = [round(x, 2) for x in np.arange(0.0,1.01,0.01)] 

allowed_weight_classes = ["unweighted","linear","PHRED","One_minus_proportion","One_over_proportion","One_over_mutation_rate"]

o_over_e_col = ["Unweighted_O_over_E","Linear_Weighted_O_over_E","PHRED_Weighted_O_over_E","One_minus_prop_Weighted_O_over_E","One_over_prop_Weighted_O_over_E","One_over_mr_Weighted_O_over_E"]

#---------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
#---------------------------------------------------------------------------------------------------------------------------------

def add_constraint_scores(sub_p):

    p = sub_p.add_parser("score",
                         help = "Calculate the O/E and Percentile constraint scores",
                         formatter_class=argparse.RawDescriptionHelpFormatter,
                         description=("\n\t*****************************************\n"
                                      "\t* O/E and Percentile Constraint Scoring *\n"
                                      "\t*****************************************\n\n"
                                      "\tCalculate the Observed over Expected (O/E) score and\n"
                                      "\tthe Percentile constraint scores for each constraint region\n\n"
                                      "\t - The percentile score is calculated after all regions have an O/E score\n"
                                      "\t   and represent the genome-wide constraint percentile of one region to all\n"
                                      "\t   other regions."
                         )
    )

    req = p.add_argument_group("Required Arguments")
    req.add_argument(
        "--o-and-e-scores",
        metavar = "Observed and Expected Score File", 
        required=True, 
        help="(Required) The path the observed over expected scores file."
    )

    req.add_argument(
        "--mutation-table",
        metavar="Mutation Table",
        required = True,
        help = "(Required) The mutation table used to calculate the mutation rate. The different weights used to calcluate the O/E scores will be calculated from this mutation frequency table"
    )

    req.add_argument(
        "--out-file",
        metavar="Ouput File",
        required = True,
        help = "(Required) The path and/or the name of the output file to create. This output file will containe the same content as the original file with additional columns for O/E scores and Percentile Scores"
    )

    p.add_argument(
        "--pct-rec-rate",
        metavar = "Precent Recovery Rate",
        default = 0.8,
        choices = recovery_choices,
        help = "The percent/fraction of bases of the total bases of a region that are rquired to be recovered (bases that were used for the region score). The O/E score along with the percentile score will be calculated for any region with a percent recovery rate >= this value. (Default = 0.8, meaing 80%% Recovery Rate)"
    )

    p.add_argument(
        "--remove-duplicate",
        action="store_true",
        help="Whether or not to remove duplicate gene entries if they exists. Default is set to False. This argument should not be set if there are multiple regions with scores for a single gene. This argument should be set if each region is a single gene and where multiple scores ofr a single gene is bad"
    )

    p.add_argument(
        "--pct-col-name",
        metavar="Percentile Column Name",
        default = "sc_percentile",
        help = "The name of the column to create the represents the percentile score. (Default = 'sc_percentile' for splicing constarint percentile)"
    )

    p.add_argument(
        "--sort-by-pos",
        action="store_true",
        help="If the `--sort-by-pos` argument is set, the output file will be sorted by the chromosome and genomic positions. If this argument is not set then the output will be sorted by increasing percentile score. (Default = sort by increasing percentile score)"
    )

    p.set_defaults(func=constraint_scores)

#---------------------------------------------------------------------------------------------------------------------------------
## Functions/Methods
#---------------------------------------------------------------------------------------------------------------------------------

def convert_scores_to_percentiles(query_df, score_column, percentile_column_name, invert_percentiles = False):
    """
    convert_scores_to_percentiles
    =============================
    Method used to convert a range of scores to a range of percentiles from 0.0 to 1.0. Most often, the higher 
     the percentile the better the score. This method will use the range of scores from a "scores" column in a 
     dataframe and based on the range provide a percentile of which each score fits into the range. (From 0 to 1)
    
    Paremeters:
    -----------
    1) query_df:               (pandas DF)  A dataframe with a "scores" column in it. 
    2) score_column:           (str) The column in the df that represents the "scores" that will be used to generate percentiles.
    3) percentile_column_name: (str) The name of the percentile column to create
    4) invert_percentiles:     (bool) Whether or not to invert the percentiles. That is, percentiles are generated 
                                      with the largest scores getting the highest percentiles and the lowest scores
                                      getting the lowest percentiles. If your better scores are smaller and your worse 
                                      scores are higher, you can invert the percentiles so that the low scores get the 
                                      high percentiles and the high scores get the low percentiles.
                                      (Default = False)
    
    Returns:
    ++++++++
    1) (Pandas DataFrame) The query df with the 'percentile_column_name' column with percentile scores added
    """
    
    n_scores = query_df.shape[0]

    percentile_scores = []

    ## Sort query df by scores
    ## Scores are sorted based on whether or not the scores need to be invereted 
    sorted_query_df = query_df.sort_values(by=[score_column], ascending=False) if invert_percentiles else query_df.sort_values(by=[score_column])
    sorted_query_df.reset_index()

    ## iterate over the size of the df and get a percentile score for each row based on the index. 
    ## The index represents the score sorted based on the score value and whether or not to invert the percentile
    ## This index represents the percentile score the specific index is within the sorted distrubtion 
    ## The index divided by the total number of score values gives you the percentile of that specific index
    for i in range(1,n_scores + 1):
        
        percentile_scores.append(float(i)/float(n_scores))

    sorted_query_df[percentile_column_name] = percentile_scores

    return(sorted_query_df)


def old_convert_scores_to_percentiles(query_df, score_column, invert_percentiles = False):
    """
    old_convert_scores_to_percentiles
    =============================
    NOTE: This method is much slower 

    Method used to convert a range of scores to a range of percentiles from 0.0 to 1.0. Most often, the higher 
     the percentile the better the score. This method will use the range of scores from a "scores" column in a 
     dataframe and based on the range provide a percentile of which each score fits into the range. (From 0 to 1)
    
    Paremeters:
    -----------
    1) query_df:           (pandas DF)  A dataframe with a "scores" column in it. 
    2) score_column:       (str)  The column in the df that represents the "scores" that will be used to generate percentiles.
    3) invert_percentiles: (bool) Whether or not to invert the percentiles. That is, percentiles are generated 
                                  with the largest scores getting the highest percentiles and the lowest scores
                                  getting the lowest percentiles. If your better scores are smaller and your worse 
                                  scores are higher, you can invert the percentiles so that the low scores get the 
                                  high percentiles and the high scores get the low percentiles.
                                  (Default = False)
    
    Returns:
    ++++++++
    1) (list) A list of percentiles, where each percentile corresponds to a score in the query_df. The index 
               of each item in the list corresponds to the row in the query_df of the score the percentile 
               was generated from. 
    """
    
    ## set the multiplier 
    multiplier = -1 if invert_percentiles else 1
    
    ## Get a list of the all scores, where each score is multiplied by the multiplier 
    scores_list = (query_df[score_column] * multiplier).to_list()
    
    ## Convert the scores in the score_column to percentiles. 
    percentile_list = query_df[score_column].apply(lambda x: (stats.percentileofscore(scores_list,(x * multiplier))/100)).to_list()
    
    return(percentile_list)


def calculate_o_over_e(lined, observed_column_suffix, expected_column_suffix, one_exon_mean = 60, weight_dict = {}, col_list = []):
    """
    calculate_o_over_e
    ===================
    Method to calculate the Observed over Expected (O/E) score for a specific region, a row in a pandas data frame, based on 
     the observed scores and expected scores split across different reference allele delta bins. 

    NOTE: This method is set up to be used as the function input when running a pandas data frame apply function.

    Detal bins: bins that represent different delta scores by which the O and E scores are seperated/binned by. (Splice AI Delta bins)

    Reference Allele: Each detal bin is further seperated by the reference allele. That is, for each delta bin category there are four 
     bin scores based on the a reference allele of A, C, G, and T.


    O/E Equation:

     SUM( weight[x] * (( Observed Counts[i] - Expected Counts[i] ) / ( Expected Counts[i] )) )

        where x = one of the possible weights. For example, the linear weight uses the top score for a delta score bin.

        where i = a reference speicifc delta score bin for all reference score bins. (If 6 delta bins, 6 bins * 4 ref alleles = 24 ref allele delta bins) 

    Parameters:
    -----------
    1) df_row:    (Pandas DataFrame)  A row from a pandas data frame. (Usually from a pandas .apply function)
    2) observed_column_suffix: (str)  The suffix of the column that represents the observed counts in the pandas data frame
    3) expected_column_suffix: (str)  The suffix of the column that represents the expected counts in the pandas data frame
    4) one_exon_mean:          (int)  The mean value to set for the distribution of regions with a single exon. (Default = 60)
    5) weight_class:           (str)  The weigthing class to use. Must be on of: "unweighted","linear","PHRED","One_minus_proportion","One_over_proportion","One_over_mutation_rate"  
    7) weight_dict:           (dict)  A 2d dictionary with keys as a weight class, values as the Weights for that class based on delta bin and ref.

    Returns:
    ++++++++
    1) (float) The Observed / Expected (O/E) score for the current region (pandas data frame row)
    """

    ## Check the weight class
#    if weight_class not in allowed_weight_classes:
#        print("!!ERROR!! The '{}' weighting class is not allowed. Please use one of: \n\t{}".format(", ".join(allowed_weight_class)))
#        sys.exit(1)

    ## Get the Observed and Expected counts for each reference allele specific detal bin
    o_over_e_scores = [[] for _ in range(len(allowed_weight_classes))]
    for delta_bin in by_ref_delta_score_bins:
        expected = float(lined["{}_{}".format(delta_bin, expected_column_suffix)])
        observed = float(lined["{}_{}".format(delta_bin, observed_column_suffix)])

        for weight_index, weight_class in enumerate(allowed_weight_classes):

            ## Identify the weight to apply
            weight = (1 if weight_class == "unweighted" else 
                        float(delta_bin.strip().split("-")[1].strip().split("_")[0]) if weight_class == "linear" else
                        weight_dict[weight_class][delta_bin])

            ## Multiply O-E/E by weight
            o_over_e_scores[weight_index].append((weight * (((observed - expected)) / (expected))))

    ## Any gene that has a max number of exons <= 1 should be considered unconstrained 
    ### These genes are artisally set to a mean of 30, with their O/E score distributed around that mean
    #### This allows the genes O/E score to be used, but all of these genes will be the most unconstarined genes

    for i,score in enumerate(o_over_e_scores):

        lined[col_list[i]] = (sum(score)  + one_exon_mean) if float(lined["max_exon_number"]) <= 1 else sum(score) 


def get_weights(mutation_table):
    """
    get_weights
    =================
    Method to calcuate different weights to add to calculate the O/E score with using a mutation frequency table. This weight can 
     be used as the O/E scaling factor.
    
     1) The PHRED scores is based on the marginal proportions of sites 
         in the mutation frequency table. That is, for each reference allele, was is the proportion of sites
         with that reference allele and SpliceAI score bin. The PHRED equation (-10 * log10(Proportion)) is 
         used on the marginal proportions to come up with the PHRED-like weight. 
    2) One minus proportion (1 - proprotion) provides a weighted based on the proportion of sites at a given
        reference allele and spliceAI bin
    3) One over proportion (1/proportion) provides a scaled weighted based on the proportion of sites
    4) One over mutation rate (1/mutation rate) provides a scaled weighted based on the mutation rate 


    Parameters:
    -----------
    1) mutation_table: (str) The file path to the mutation table to use to calcualte that PHRED weight 

    Returns:
    ++++++++
    1) (dict) A dictionary with keys as weight class, value as a second dictionary with keys as {delta score bin}_reference alleles,m and values as the weight
                weight classes:
                    PHRED
                    One_minus_proportion
                    One_over_proportion
                    One_over_mutation_rate
    """

    from math import log10
    
    header_index = 0
    try: 
        with io.open(mutation_table, "rt", encoding = "utf-8") as mrt: 
            for i,line in enumerate(mrt):
                if line[0] == "#": 
                    header_index = i
                elif line[0] != "#": 
                    break
    except IOError as e:
        print("\n!!ERROR!! Unable to read the mutation rate table")
        print(str(e))
        sys.exit(1)

    ## Load the nss table into a pandas df
    mr_table = pd.read_csv(mutation_table, sep="\t", index_col=False, header = header_index)
    mr_table = mr_table.rename(columns = {"#ref":"ref"})

    ## Get the sum of counts for each reference allele and SpliceAI score bin combination
    by_ref_counts = mr_table.groupby(["delta_score_bin","ref"]).agg({"zerotons":sum, "non_zerotons":sum}).reset_index()
    by_ref_counts["total"] = by_ref_counts.zerotons + by_ref_counts.non_zerotons

    ## Create a dictionary of total marginal counts by reference allele
    total_dict = by_ref_counts.groupby("ref").total.sum().reset_index().set_index("ref").T.to_dict() 

    ## Add the by reference marginal proportions for each reference allele and SpliceAI score combination 
    by_ref_counts["ref_marginal_proportion"] = by_ref_counts.apply(lambda x: x.total/total_dict[x.ref]["total"], axis = 1)

    ## Calculate the PHRED-like weight using the marginal proportions 
    ### -10 * log10(Proportion)
    by_ref_counts["PHRED_Weight"] = by_ref_counts.apply(lambda x: (-10 * (log10(x.ref_marginal_proportion))), axis =1)

    weights_dict = defaultdict(lambda: (defaultdict(float)))

    print("\n\tPHRED WEIGHTS:")
    print("\t==============")
    print("\n\t SpliceAI_Ref\tPHRED Weight")
    print("\t ------------\t------------")
    for row in by_ref_counts.itertuples():
        
        print("\t {}_{}:\t{}".format(row.delta_score_bin, row.ref, row.PHRED_Weight))
        weights_dict["PHRED"]["{}_{}".format(row.delta_score_bin, row.ref)] = row.PHRED_Weight


    ## ! - marginal proportion 
    by_ref_counts["One_min_proportion"] = 1 - by_ref_counts.ref_marginal_proportion

    print("\n\t1 - proportion Weights:")
    print("\t=======================")
    print("\n\t SpliceAI_Ref\t1 - Proportion")
    print("\t ------------\t--------------")
    for row in by_ref_counts.itertuples():
        
        print("\t {}_{}:\t{}".format(row.delta_score_bin, row.ref, row.One_min_proportion))
        weights_dict["One_minus_proportion"]["{}_{}".format(row.delta_score_bin, row.ref)] = row.One_min_proportion


    ## 1 / marginal proportion 
    by_ref_counts["One_over_proportion"] = 1/by_ref_counts.ref_marginal_proportion

    print("\n\t1 / proportion Weights:")
    print("\t=======================")
    print("\n\t SpliceAI_Ref\t1 / Proportion")
    print("\t ------------\t--------------")
    for row in by_ref_counts.itertuples():
        
        print("\t {}_{}:\t{}".format(row.delta_score_bin, row.ref, row.One_over_proportion))
        weights_dict["One_over_proportion"]["{}_{}".format(row.delta_score_bin, row.ref)] = row.One_over_proportion

    ## 1 / mutation rate
    by_ref_counts["mutation_rate"] = by_ref_counts.non_zerotons / by_ref_counts.total
    by_ref_counts["One_over_mutation_rate"] = 1 / by_ref_counts.mutation_rate

    print("\n\t1 / mutation rate Weights:")
    print("\t==========================")
    print("\n\t SpliceAI_Ref\t1 / MutationRate")
    print("\t ------------\t----------------")
    for row in by_ref_counts.itertuples():
        
        print("\t {}_{}:\t{}".format(row.delta_score_bin, row.ref, row.One_over_mutation_rate))
        weights_dict["One_over_mutation_rate"]["{}_{}".format(row.delta_score_bin, row.ref)] = row.One_over_mutation_rate


    return(weights_dict)


#---------------------------------------------------------------------------------------------------------------------------------
## Main
#---------------------------------------------------------------------------------------------------------------------------------

def constraint_scores(parser, args):

    print("\n\n\t*****************************************")
    print("\t* O/E and Percentile Constraint Scoring *")
    print("\t****************************************\n\n")


    print(("\nInput Arguments:"
           "\n================"
           "\n - o-and-e-scores:           {}"
           "\n - mutation-table            {}" 
           "\n - out-file:                 {}"
           "\n - pct-rec-rate:             {}"
           "\n - remove-duplicate:         {}"
           "\n - pct-col-name:             {}"
           "\n - sort-by-pos:              {}"

           "\n"
           ).format(args.o_and_e_scores, 
                    args.mutation_table,
                    args.out_file,
                    args.pct_rec_rate,
                    args.remove_duplicate,
                    args.pct_col_name,
                    "Output will be sorted by genomic positions" if args.sort_by_pos else "Output will be sorted by increasing percentile",
                    )
    )

    ## Get different weights 
    print("\nCalculating  weights to use for O/E scoring")
    print("\n  NOTE: Only calculated from the zeroton model")
    
    weights_dict = get_weights(args.mutation_table)

    ## load the O and E Score file into a pandas data frame
    print("\nReading data from: {}".format(args.o_and_e_scores))

    try:
        fh = io.open(args.o_and_e_scores, "rt", encoding="utf-8")
    except IOError as e:
        print("\n\n!!ERROR!! Unable to read '{}'. Please correct the error and try again.".format(args.o_and_e_scores))
        print(str(e))
        sys.exit(1)

    
    header = []
    region_count = 0
    filtered_region_count = 0
    score_dict = dict()
    good_keys = set()

    print("\nParsing O/E scores and applying filters")
    for line in fh:
        
        if line[0] == "#":
            header = line.strip().replace("#","").split("\t")
            continue

        line_dict = dict(zip(header,line.strip().split("\t")))
        region_count += 1

        ## Filters 
        ## 1. Remove any regions with no expecatation scores. (These regions did not received o or e scores)
        if float(line_dict["zeroton_expectation_sum"]) <= 0.0:
            continue

        ## 2. Identify the percent/fraction of recovered bases per region
        total_positions = "total_positions" if "total_positions" in line_dict else "total_region_positions" if "total_region_positions" in line_dict else "total_gene_positions"
        line_dict["fraction_recovered"] = float(line_dict["positions_considered"]) / float(line_dict[total_positions])

        ## 3. Keep only positions that are >= the percent/fraction reovery rate
        if line_dict["fraction_recovered"] < args.pct_rec_rate:
            continue

        filtered_region_count += 1

        calculate_o_over_e(lined = line_dict, 
                           observed_column_suffix = "zeroton_observed", 
                           expected_column_suffix = "zeroton_expected", 
                           one_exon_mean = 2200,
                           weight_dict = weights_dict, 
                           col_list = o_over_e_col)

        ## Add scores to score_dict. Only include the score columns
        dict_key = "{}:{}-{}:{}".format(line_dict["chrom"],
                                        line_dict["region_start"],
                                        line_dict["region_end"],
                                        line_dict["gene_id"]) 

        score_dict[dict_key] = {key:value for key,value in line_dict.items() if key in o_over_e_col}

        good_keys.add(dict_key)

    fh.close()

    print("\n\tNumber of regions before filtering: {}".format(region_count))
    print("\n\tNumber of regions after filtering: {}".format(filtered_region_count))


    ## Convert score dict into pandas DF
    o_and_e_df = pd.DataFrame.from_dict(score_dict, orient="index")
    
    ## Re-order columns
    o_and_e_df = o_and_e_df[o_over_e_col]

    del score_dict


    ## Calculate Percentile Score
    print("\nTransforming O/E scores into percentiles")

    print("\n\tunweighted")
    o_and_e_df = convert_scores_to_percentiles(query_df = o_and_e_df, 
                                               score_column = "Unweighted_O_over_E", 
                                               percentile_column_name = args.pct_col_name,
                                               invert_percentiles = True)

    print("\n\tlinear")
    o_and_e_df = convert_scores_to_percentiles(query_df = o_and_e_df, 
                                               score_column = "Linear_Weighted_O_over_E", 
                                               percentile_column_name = "weighted_%s" %args.pct_col_name,
                                               invert_percentiles = True)

    print("\n\tPHRED")
    o_and_e_df = convert_scores_to_percentiles(query_df = o_and_e_df, 
                                               score_column = "PHRED_Weighted_O_over_E", 
                                               percentile_column_name = "PHRED_weighted_%s" %args.pct_col_name,
                                               invert_percentiles = True)

    print("\n\t1 - Proportion")
    o_and_e_df = convert_scores_to_percentiles(query_df = o_and_e_df, 
                                               score_column = "One_minus_prop_Weighted_O_over_E", 
                                               percentile_column_name = "one_minus_prop_weighted_%s" %args.pct_col_name,
                                               invert_percentiles = True)

    print("\n\t1 / Proportion")
    o_and_e_df = convert_scores_to_percentiles(query_df = o_and_e_df, 
                                               score_column = "One_over_prop_Weighted_O_over_E", 
                                               percentile_column_name = "one_over_prop_weighted_%s" %args.pct_col_name,
                                               invert_percentiles = True)

    print("\n\t1 / Mutation Rate")
    o_and_e_df = convert_scores_to_percentiles(query_df = o_and_e_df, 
                                               score_column = "One_over_mr_Weighted_O_over_E", 
                                               percentile_column_name = "one_over_mr_weighted_%s" %args.pct_col_name,
                                               invert_percentiles = True)


    pctl_cols = [args.pct_col_name, 
                "weighted_%s" %args.pct_col_name, 
                "PHRED_weighted_%s" %args.pct_col_name, 
                "one_minus_prop_weighted_%s" %args.pct_col_name,
                "one_over_prop_weighted_%s" %args.pct_col_name,
                "one_over_mr_weighted_%s" %args.pct_col_name]


    ## Convert all values to strings
    o_and_e_df = o_and_e_df.astype(str)

    o_and_e_df_dict = o_and_e_df.T.to_dict("list")

    print("\nCreating output file '{}'".format(args.out_file))



    try:
        fh = io.open(args.o_and_e_scores, "rt", encoding="utf-8")
    except IOError as e:
        print("\n\n!!ERROR!! Unable to read '{}'. Please correct the error and try again.".format(args.o_and_e_scores))
        print(str(e))
        sys.exit(1)

    
    header = []

    with open(args.out_file, "w") as out:
        
        for line in fh:
            
            if line[0] == "#":
                header = line.strip().replace("#","").split("\t")
                out.write("#" + "\t".join(header + o_over_e_col + pctl_cols) + "\n")  
                continue

            line_list = line.strip().split("\t")

            line_dict = dict(zip(header,line_list))

            dict_key = "{}:{}-{}:{}".format(line_dict["chrom"],
                                            line_dict["region_start"],
                                            line_dict["region_end"],
                                            line_dict["gene_id"]) 
            
            ## Skip bad keys
            if dict_key not in good_keys:
                continue


            ## Write line out
            out.write("\t".join(line_list + o_and_e_df_dict[dict_key]) + "\n")

        fh.close()


    print("\nDONE\n")
    
