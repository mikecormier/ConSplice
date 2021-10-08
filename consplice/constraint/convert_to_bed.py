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

def add_to_bed(sub_p):
    
    p = sub_p.add_parser("to-bed",
                         help = "Convert the 1-based scored ConSplice txt file to a 0-based bed file",
                         formatter_class=argparse.RawDescriptionHelpFormatter,
                         description=("\n\t**********************\n"
                                      "\t* ConSplice - to bed *\n"
                                      "\t**********************\n\n"
                                      "\tConvert the scored ConSplice txt file, where the genomic positions are\n"
                                      "\t1-based like a vcf file, to a bed file, where the genomic positions will\n"
                                      "\tbe 0-based"
                        )

    )


    req = p.add_argument_group("Required Arguments")
    req.add_argument(
        "--score-file",
        metavar = "Regional Score File", 
        required=True, 
        help="(Required) The path to the scored ConSplice txt file to convert to 0-based bed file"
    )

    req.add_argument(
        "--out-file",
        metavar = "Output file",
        required = True,
        help = "(Required) The name of the output bed file to create."
    )

    p.add_argument(
        "--out-type",
        metavar="Output file type",
        choices = ["bed","bedgz"],
        default = "bed",
        help="(Optional) The output file type. Whether the output file should be a normal bed file or a bgzipped bed file. Choices = 'bed' or 'bedgz'. Default = 'bed'"
    )

    p.add_argument(
        "--sort",
        action = "store_false",
        help="(Optional) Whether to sort the output file or not. Deafult = True."
    )


    p.set_defaults(func=to_bed)


#---------------------------------------------------------------------------------------------------------------------------------
## Main
#---------------------------------------------------------------------------------------------------------------------------------

def to_bed(parser, args):


    print("\n\t**********************")
    print("\t* ConSplice - to bed *")
    print("\t**********************\n\n")


    print(("\nInput Arguments:"
           "\n================"
           "\n - score-file:           {}"
           "\n - out-file:             {}"
           "\n - out-type:             {}"
           "\n - sort:                 {}"
           ).format(args.score_file, 
                    args.out_file, 
                    args.out_type, 
                    args.sort,
                    )
    )






