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


#---------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
#---------------------------------------------------------------------------------------------------------------------------------

def add_select_score(sub_p):
    
    p = sub_p.add_parser("select-score",
                         help = "Select an O/E and matching Percentile score to filter on and remove all other non-essential columns after ConSplice scoring",
                         formatter_class=argparse.RawDescriptionHelpFormatter,
                         description=("\n\t**************************\n"
                                      "\t* Select ConSplice Score *\n"
                                      "\t**************************\n\n"
                                      "\tReduce unnecessary columns in the scored ConSplice file by filter\n"
                                      "\tout non-essential columns while keeping the desired O/E and Percentile\n"
                                      "\tConSplice scores.\n"
                                      "\tAny extra O/E scores and Percentile scores will be removed along with\n"
                                      "\tall other non-essential columns in the output file."
                        )

    )


    req = p.add_argument_group("Required Arguments")
    req.add_argument(
        "--score-file",
        metavar = "Regional Score File", 
        required=True, 
        help="(Required) The path to the scored ConSplice file to filter"
    )

    req.add_argument(
        "--o-over-e-field",
        metavar = "O/E Field",
        required = True,
        help = "(Required) The name of the O/E score field/column within the scored ConSplice file to select/keep."
    )

    req.add_argument(
        "--pctl-field",
        metavar = "Percentile Field",
        required = True,
        help = "(Required) The name of the Percentile score field/column within the scored ConSplice file to select/keep."
    )

    req.add_argument(
        "--out-file",
        metavar = "Output file",
        required = True,
        help = "(Required) The name of the output file to create."
    )

    p.add_argument(
        "--new-oe-name",
        metavar="New O/E Column Name",
        default = "ConSplice_O/E",
        help="(Optional) The new name of the O/E score field/column in the output file. Default = 'ConSplice_O/E'"
    )

    p.add_argument(
        "--new-pctl-name",
        metavar="New Percentile Score Column Name",
        default = "ConSplice_Percentile",
        help="(Optional) The new name of the Percentile score field/column in the output file. Default = 'ConSplice_Percentile'"
    )


    p.set_defaults(func=select_scores)


#---------------------------------------------------------------------------------------------------------------------------------
## Main
#---------------------------------------------------------------------------------------------------------------------------------

def select_scores(parser, args):

    print("\n\nUPDATE SCRIPT TO ADD A HEADER LINE WITH WHICH O/E and PCTL Score WAS KEPT\n\n")


    print("\n\nUPDATE SCRIPT TO READ and WRITE In Place rather than loading in with pandas\n\n")


    print("\n\nUPDATE SCRIPT TO select scores for both gene and region files. (Currently only works for region files)\n\n")


    print("\n\t**************************")
    print("\t* Select ConSplice Score *")
    print("\t**************************\n\n")



    print(("\nInput Arguments:"
           "\n================"
           "\n - score-file:           {}"
           "\n - o-over-e-field:       {}"
           "\n - pctl-score:           {}"
           "\n - out-file:             {}"
           "\n - new-o-over-e-name:    {}"
           "\n - new-pctl-name:        {}"
           ).format(args.score_file, 
                    args.o_over_e_field, 
                    args.pctl_field, 
                    args.out_file, 
                    args.new_o_over_e_name, 
                    args.new_pctl_name
                    )
    )

    print("\n> Reading score file")
    ## Read scores file into dataframe
    scores_df = pd.read_csv(args.score_file, sep = "\t", index_col = False)

    print("\n> Filtering columns.")
    
    print("\n\tKeeping the '{}' O/E and '{}' Percentile score fields".format(args.o_over_e_field, args.pctl_field))

    print("\n\tRenaming the '{}' O/E score field to '{}".format(args.o_over_e_field, args.new_o_over_e_name))

    print("\n\tRenaming the '{}' Percentile score field to '{}'".format(args.pctl_field, args.new_pctl_name))

    scores_df = scores_df.rename(columns = {args.o_over_e_field: args.new_o_over_e_name, args.pctl_field: args.new_pctl_name})

    scores_df["region_size"] = scores_df.region_end - (scores_df.region_start - 1) ## offset by 1 to count for the first position and last position, not just the difference

    scores_df = scores_df[["#chrom",
                           "region_start",
                           "region_end",
                           "txStart",
                           "txEnd",
                           "strand",
                           "transcript_id",
                           "gene_id",
                           "gene_symbol",
                           "region_size",
                           args.new_o_over_e_name,
                           args.new_pctl_name]]

    print("\n> Sorting region by alphanumeric chromosome and start position")
    scores_df = scores_df.sort_values(by = ["#chrom", "region_start"])

    print("\n> Wriring output to '{}'".format(args.out_file))
    scores_df.to_csv(args.out_file, sep="\t", index=False)

    print("\n> DONE")


if __name__ == "__main__":
    sys.exit(main() or 0)



