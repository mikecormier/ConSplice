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



#---------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
#---------------------------------------------------------------------------------------------------------------------------------

def add_ml_scoring(sub_p):

    p = sub_p.add_parser("score",
                         help = "Score variants using ConSpliceML",
                         formatter_class=argparse.RawDescriptionHelpFormatter,
                         description = ("\n\t*********************\n\t* ConSpliceML Score *\n\t*********************\n\n"
                                        "\tScore variants using the trained ConSpliceML model"
                                       ),
    )


    req = p.add_argument_group("Required Arguments")
    req.add_argument(
        "--input-type",
        metavar = "Input File Type", 
        choices = ["vcf","txt","bed"],
        required=True, 
        help="(Required) The intput file type of the variants to score. ConSpliceML can score either a vcf file or a txt/bed file. The output file type will match the input file type where the score will be added to the INFO field of a vcf file and the score will be added as a column to the txt/bed file. NOTE: vcf and txt files are assumed to be 1-based while bed files are assumed to be 0-based"
    )

    req.add_argument(
        "--input-file",
        metavar = "Input File",
        required = True,
        help = "(Required) The input variant file to score"
    )

    req.add_argument(
        "--cs-model",
        metavar = "ConSplice Model",
        required = True,
        help = "(Required) The path to the directory that contains the ConSpliceML Model."
    )

    p.add_argument(
        "--add-consplice",
        action = "store_true",
        help = "(Optional) Whether or not to add the regional constrained splicing score to the output file. Default = False"
    )

    p.set_defaults(func=rf_score)


#---------------------------------------------------------------------------------------------------------------------------------
## Functions/Methods
#---------------------------------------------------------------------------------------------------------------------------------




#---------------------------------------------------------------------------------------------------------------------------------
## Main
#---------------------------------------------------------------------------------------------------------------------------------

def rf_score(parser, args):


    print("\n\n\t*************************")
    print("\t** ConSpliceML - Score **")
    print("\t*************************\n\n")

    print(("\nInput Arguments:"
           "\n================"
           "\n - intput-type:           {}"
           "\n - intput-file:           {}"
           "\n - cs-model:              {}"
           "\n - add-consplice:         {}"
           "\n"
           ).format(args.input_type, 
                    args.input_file, 
                    args.cs_model, 
                    args.add_consplice,
                    )
    )



