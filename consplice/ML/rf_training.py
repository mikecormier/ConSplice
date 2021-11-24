from __future__ import print_function
import sys
import os
import io
from pathlib import Path
import argparse 
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
import joblib
import yaml
import datetime
from .ml_utils import ConSpliceML_train, patho_label, trained_model, ml_info_file


#---------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
#---------------------------------------------------------------------------------------------------------------------------------

def add_training(sub_p):

    p = sub_p.add_parser("train",
                         help = "Train a Random Forest model using ConSplice.",
                         formatter_class=argparse.RawDescriptionHelpFormatter,
                         description = ("\n\t****************************\n"
                                        "\t** ConSpliceML - Training **\n"
                                        "\t****************************\n\n"
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
        help="(Required) The path and file name to the pathogenic variants for training. ConSplice, SpliceAI, and SQUIRL scores should be in this file. Requires a tab-delimited file. NOTE: This file should be split into the training set prior to using it here."
    )

    req.add_argument(
        "--benign-set",
        metavar = "Benign Splicing set", 
        required=True, 
        help="(Required) The path and file name to the benign variants for training. ConSplice, SpliceAI, and SQUIRL scores should be in this file. Requires a tab-delimited file. NOTE: This file should be split into the training set prior to using it here."
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
        "--output-dir",
        metavar = "Output Directory",
        default = "ConSpliceML_Model",
        help = "(Optional) The path and name of the directory to crate and store the final trained model in. The default output dir will be created in the current working directory under the name 'ConSpliceML_Model'"
    )

    p.add_argument(
        "--random-state",
        metavar = "Random State",
        default = 156498,
        type = int,
        help = "(Optional) The random state to use when training the model. (Default = 156498)"
    )

    p.add_argument(
        "--n-dt",
        metavar = "N decision trees",
        default = 1000,
        type = int,
        help = "(Optional) The number of decision trees to use when creating the random forest  (Default = 1000)"
    )

    p.set_defaults(func=rf_training)



#---------------------------------------------------------------------------------------------------------------------------------
## Main
#---------------------------------------------------------------------------------------------------------------------------------

def rf_training(parser, args):


    print("\n\n\t****************************")
    print("\t** ConSpliceML - Training **")
    print("\t****************************\n\n")


    ## Check random state
    if not isinstance(args.random_state, int):
        print("\n!!ERROR!! The random state needs to be a positive integer")
        sys.exit(1)

    if args.random_state < 0:
        print("\n!!ERROR!! The random state needs to be a positive integer")
        sys.exit(1)

    ## Check n decision tress
    if not isinstance(args.n_dt, int):
        print("\n!!ERROR!! The n-dt parameter requires a positive integer")
        sys.exit(1)

    if args.n_dt < 1:
        print("\n!!ERROR!! The n-dt parameter requires a positive integer")
        sys.exit(1)


    ## Check that the patho and benign files exist
    if not os.path.exists(os.path.abspath(args.patho_set)) or not os.path.isfile(os.path.abspath(args.patho_set)):
        print("\n!!ERROR!! The pathogenic training set file does not exists or is not a file. Please correct the error and try again")
        sys.exit(1)

    if not os.path.exists(os.path.abspath(args.benign_set)) or not os.path.isfile(os.path.abspath(args.benign_set)):
        print("\n!!ERROR!! The benign training set file does not exists or is not a file. Please correct the error and try again")
        sys.exit(1)

    output_dir = os.path.abspath(args.output_dir)

    print(("\nInput Arguments:"
           "\n================"
           "\n - patho-set:            {}"
           "\n - benign-set:           {}"
           "\n - consplice-col:        {}"
           "\n - spliceai-col:         {}"
           "\n - squirls-col:          {}"
           "\n - output-dir:           {}"
           "\n - random-state:         {}"
           "\n - n-dt:                 {}"
           "\n"
           ).format(args.patho_set, 
                    args.benign_set, 
                    args.consplice_col, 
                    args.spliceai_col, 
                    args.squirls_col,
                    output_dir,
                    args.random_state,
                    args.n_dt,
                    )
    )

    ## Load patho set
    print("\nLoading pathogenic training set")
    patho_df = pd.read_csv(args.patho_set, sep = "\t") 
    patho_df[patho_label] = 1

    print("\n\tloaded {} training pathogenic samples.".format(patho_df.shape[0]))

    ## Load benign set
    print("\nLoading bening training set")
    benign_df = pd.read_csv(args.benign_set, sep = "\t") 
    benign_df[patho_label] = 0

    print("\n\tloaded {} training benign samples.".format(benign_df.shape[0]))

    ## gather features
    print("\nExtracting features")

    print("\n\tExtracting features using '{}', '{}', and '{}' feature names".format(args.consplice_col, args.spliceai_col, args.squirls_col))
    combined_df = pd.concat([patho_df[[args.consplice_col, args.spliceai_col, args.squirls_col, patho_label]], 
                             benign_df[[args.consplice_col, args.spliceai_col, args.squirls_col, patho_label]]])

    print("\n\tRemoving features with missing values")
    combined_df = combined_df.loc[~(combined_df[args.consplice_col].isna()) & 
                                  ~(combined_df[args.spliceai_col].isna()) & 
                                  ~(combined_df[args.squirls_col].isna()) ] 

    print("\n\tshuffling samples")
    combined_df = combined_df.sample(frac=1, random_state = args.random_state)

    print("\nTraining the Random Forest using: \n\n\t- {} decision trees \n\n\t- {} pathogenic samples \n\n\t- {} benign samples".format(args.n_dt, combined_df.loc[combined_df[patho_label] == 1].shape[0], combined_df.loc[combined_df[patho_label] == 0].shape[0]))
    rf = ConSpliceML_train(training_df = combined_df,
                           feature_col_names = [args.consplice_col, args.spliceai_col, args.squirls_col],
                           label_col = patho_label,
                           n_estimators = args.n_dt,
                           random_state = args.random_state)

    print("\nPreparing to save")
    ## check if the output directory exists
    ### Create it if it does not
    if not os.path.exists(os.path.abspath(output_dir)):
        print("\n\tOutput directiory '{}' does not exists. Creating output directory".format(output_dir))
        Path(output_dir).mkdir(parents = True)
        if os.path.exists(output_dir) and os.path.isdir(output_dir):
            print("\n\tSuccessfully created the output direcotry '{}'".format(output_dir))
        else:
            print("\n!!ERROR!! Unabel to create the output directory")
            sys.exit(1)

    else:
        if not os.path.isdir(os.path.abspath(output_dir)):
            print("\n!!ERROR!! The '{}' name exists but is not a directory. Please change the directory path/name or fix the problem".format(output_dir))
            sys.exit(1)

        else:
            print("\n\tThe output directory '{}' exists. The model will be written to this directory. Any previous model in this directory will be overwritten".format(output_dir))

    ## Save Random Forest Model 
    joblib.dump(rf, os.path.join(output_dir,trained_model))

    ## Save meta info to yaml file in the Model directory 
    info_dict = {"Features":{"ConSplice":args.consplice_col,"SpliceAI":args.spliceai_col,"SQUIRLS":args.squirls_col},
                 "Feature Order":["ConSplice","SpliceAI","SQUIRLS"],
                 "Random State":args.random_state,
                 "N Decision Trees":args.n_dt,
                 "Patho dataset": args.patho_set,
                 "Benign dataset": args.benign_set,
                 "N patho samples": combined_df.loc[combined_df[patho_label] == 1].shape[0],
                 "N benign samples": combined_df.loc[combined_df[patho_label] == 0].shape[0],
                 "Saved Model": trained_model,
                 "Saved Path": os.path.join(os.path.abspath(output_dir), trained_model),
                 "Date Created": datetime.datetime.now().strftime("%m-%d-%Y")}

    with open(os.path.join(output_dir, ml_info_file), "w", encoding = "utf-8") as fh:
        
        yaml.dump(info_dict, fh)

    print("\nThe ConSpliceML model as been saved to '{}'".format(os.path.join(output_dir, trained_model)))

    print("\nDONE\n")
