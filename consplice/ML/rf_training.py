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

def add_training(sub_p):

    p = sub_p.add_parser("train",
                         help = "Train a Random Forest model using ConSplice.",
                         formatter_class=argparse.RawDescriptionHelpFormatter,
                         description = ("\n\t************************\n\t* ConSpliceML Training *\n\t************************\n\n"
                                        "\tTrain a Random Forest using the regional ConSplice model,"
                                        "\n\tSpliceAI alternative splicing prediction scores, and"
                                        "\n\tSQUIRLS alternative splicing prediction scores."
                                       ),
    )


    req = p.add_argument_group("Required Arguments")
    req.add_argument(
        "--patho-set",
        metavar = "Pathogenic Splicing set", 
        required=True, 
        help="(Required) The path and file name to the pathogenic variants for training. ConSplice, SpliceAI, and SQUIRL scores should be in this file. NOTE: This file should be split into the training set prior to using it here."
    )

    req.add_argument(
        "--benign-set",
        metavar = "Benign Splicing set", 
        required=True, 
        help="(Required) The path and file name to the benign variants for training. ConSplice, SpliceAI, and SQUIRL scores should be in this file. NOTE: This file should be split into the training set prior to using it here."
    )

    req.add_argument(
        "--consplice-col",
        metavar = "ConSplice Column",
        required = True,
        help = "(Required) The name of the ConSplice column in the pathogenic and benign training sets"
    )

    req.add_argument(
        "--spliceai-col",
        metavar = "SpliceAI Column",
        required = True,
        help = "(Required) The name of the SpliceAI column in the pathogenic and benign training sets"
    )

    req.add_argument(
        "--squirls-col",
        metavar = "SQUIRLS Column",
        required = True,
        help = "(Required) The name of the SQUIRLS column in the pathogenic and benign training sets"
    )

    p.add_argument(
        "-output-dir",
        metavar = "Output Directory",
        default = "ConSpliceML_Model",
        help = "(Optional) The path and name of the directory to crate and store the final trained model in. The default output dir will be created in the current working directory under the name 'ConSpliceML_Model'"
    )


    p.set_defaults(func=rf_training)



#---------------------------------------------------------------------------------------------------------------------------------
## Functions/Methods
#---------------------------------------------------------------------------------------------------------------------------------




#---------------------------------------------------------------------------------------------------------------------------------
## Main
#---------------------------------------------------------------------------------------------------------------------------------

def rf_training(parser, args):


    print("\n\n\t****************************")
    print("\t** ConSpliceML - Training **")
    print("\t****************************\n\n")

    print(("\nInput Arguments:"
           "\n================"
           "\n - patho-set:            {}"
           "\n - benign-set:           {}"
           "\n - consplice-col:        {}"
           "\n - spliceai-col:         {}"
           "\n - squirls-col:          {}"
           "\n - output-dir:           {}"
           "\n"
           ).format(args.patho_set, 
                    args.benign_set, 
                    args.consplice_col, 
                    args.spliceai_col, 
                    args.squirls_col,
                    args.output_dir, 
                    )
    )



