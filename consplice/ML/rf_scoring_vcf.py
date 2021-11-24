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
import multiprocessing
from collections import defaultdict
from cyvcf2 import VCF,Writer
from .ml_utils import ConSpliceML_score, trained_model, ml_info_file,  get_alternative_gene_symbols

#---------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
#---------------------------------------------------------------------------------------------------------------------------------

def add_ml_scoring_vcf(sub_p):

    p = sub_p.add_parser("score-vcf",
                         help = "Score variants using ConSpliceML",
                         formatter_class=argparse.RawDescriptionHelpFormatter,
                         description = ("\n\t*****************************\n"
                                        "\t** ConSpliceML - score-vcf **\n"
                                        "\t*****************************\n\n"
                                        "\tScore a vcf file using the trained ConSpliceML model"
                                        "\tNOTE: The vcf file needs to be annotated with the"
                                        "\t      required splicing scores: SpliceAI, SQUIRLS, and"
                                        "\t      ConSplice. This process will fail if these annotation"
                                        "\t      are missing from the INFO field of the vcf file."
                                       ),
    )


    req = p.add_argument_group("Required Arguments")

    req.add_argument(
        "--vcf-file",
        metavar = "Input vcf File",
        required = True,
        help = "(Required) The input vcf file to score"
    )

    req.add_argument(
        "--output-file",
        metavar = "Output File",
        required = True,
        help = "(Required) The output variant file to create"
    )

    req.add_argument(
        "--alt-gene-symbol",
        metavar = "Alternative Gene Symbol File",
        default = "NA",
        help = "(Required) A file that contains mappings between a canonical gene symbol to alternative gene symbols. NOTE: The file needs to have a header!. This script is set up to use the HGNC protein-coding gene mapping file can be found in the ConSplice GitHub repo: https://github.com/mikecormier/ConSplice"  
    )

    req.add_argument(
        "--ml-model",
        metavar = "ConSpliceML Model",
        required = True,
        help = "(Required) The path to the directory that contains the ConSpliceML Model."
    )

    p.add_argument(
        "--out-type",
        metavar="Output file type",
        choices = ["vcf","vcfgz","bcf","bcfgz"],
        default = "vcf",
        help="The output file type. Whether the output file should be a normal vcf/bcf file or a bgzipped vcf/bcf file. Choices = ['vcf','vcfgz','bcf','bcfgz']. vcfgz and bcfgz are the compressed version of vcf or bcf files. Default = 'vcf'"
    )

    p.add_argument(
        "--n-cpu",
        metavar = "Number of CPUs",
        default = 3,
        type = int,
        help = "The number of CPUs to use for multi-threading during vcf parsing. (Default = 3)"
    )


    p.set_defaults(func=rf_score_vcf)


#---------------------------------------------------------------------------------------------------------------------------------
## Functions/Methods
#---------------------------------------------------------------------------------------------------------------------------------




#---------------------------------------------------------------------------------------------------------------------------------
## Main
#---------------------------------------------------------------------------------------------------------------------------------

def rf_score_vcf(parser, args):


    print("\n\n\t*****************************")
    print("\t** ConSpliceML - score-vcf **")
    print("\t*****************************\n\n")


    ## Check output file
    outfile = args.output_file
    if not args.output_file.endswith("vcf") and args.out_type == "vcf":
        outfile = args.output_file + ".vcf"

    elif not args.output_file.endswith("bcf") and args.out_type == "bcf":
        outfile = args.output_file + ".bcf"

    elif not args.output_file.endswith("vcf.gz") and args.out_type == "vcfgz":
        outfile = args.output_file + ".vcf.gz" if not args.output_file.endswith(".vcf") else args.output_file + ".gz"

    elif not args.output_file.endswith("bcf") and args.out_type == "bcfgz":
        outfile = args.output_file + ".bcf" 

    ## Check number of CPUs
    ncpus = args.n_cpu if args.n_cpu > 0 else 1
    if ncpus > multiprocessing.cpu_count(): 
        print("\n**Warning** Too many CPUs designated. Number of CPUs will be set to {}".format( multiprocessing.cpu_count() if multiprocessing.cpu_count() > 0 else 1))  
        ncpus = multiprocessing.cpu_count() if multiprocessing.cpu_count()  > 0 else 1


    print(("\nInput Arguments:"
           "\n================"
           "\n - vcf-file:              {}"
           "\n - output-file:           {}"
           "\n - alt-gene-symbol:       {}"
           "\n - ml-model:              {}"
           "\n - out-type:              {}"
           "\n - n-cpu:                 {}"
           "\n"
           ).format(args.vcf_file, 
                    outfile,
                    args.alt_gene_symbol,
                    args.ml_model, 
                    args.out_type,
                    args.n_cpu,
                    )
    )

    ## Check the ML model
    assert os.path.exists(args.ml_model), "\n!!ERROR!! The ConSpliceML model directory '{}' does not exists.\n".format(os.path.abspath(args.ml_model))

    assert os.path.isdir(args.ml_model), "\n!!ERROR!! The ConSpliceML model directory '{}' is not a directory and does not contain the ConSpliceML trained model.\n".format(os.path.abspath(args.ml_model))

    assert os.path.exists(os.path.join(os.path.abspath(args.ml_model),trained_model)), "\n!!ERROR!! The trained ConSpliceML RF model is missing from the model directory\n"

    assert os.path.exists(os.path.join(os.path.abspath(args.ml_model),ml_info_file)), "\n!!ERROR!! The metadata file for the trained ConSpliceML RF model is missing from the model directory\n"

    ## Check that the input file exists
    assert os.path.exists(os.path.abspath(args.vcf_file)), "\n!!ERROR!! The vcf file to add ConSpliceML scores to does not exists\n"

    ## Check that the alt gene symboles file
    assert os.path.exists(os.path.abspath(args.alt_gene_symbol)), "\n!!ERROR!! The alt gene symbols file does not exists\n"


    ## Check for overwritting output file
    if os.path.exists(outfile):
        print("\n\n**WARNING** The output file '{}' already exists. This file will be overwritten. To halt this process press CTRL+C\n".format(outfile)) 


    ## load the ConSpliceML model
    print("\nLoading the ConSpliceML model")
    rf = joblib.load(os.path.join(os.path.abspath(args.ml_model),trained_model))


    ## load the ConSpliceML model metadata
    print("\nLoading the ConSpliceML model metadata")
    with open(os.path.join(os.path.abspath(args.ml_model),ml_info_file)) as fh:
        
        metadata = yaml.safe_load(fh)

    feature_order = [metadata["Features"][x] for x in metadata["Feature Order"]]

    assert len(metadata) > 0, "\n!!ERROR!! The metadata file is empty or did not load properly\n"


    print("\nCreating an map of alternative gene symbols")
    ## Get a dictionary of alternative gene symbols
    alt_symbol_dict = get_alternative_gene_symbols(args.alt_gene_symbol)


    ## Open vcf file
    print("\nPreparing vcf file")
    vcf = VCF(args.vcf_file, threads=ncpus)

    ## Get the SpliceAI info and create an index
    spliceai_info_index = vcf["SpliceAI"]["Description"].strip().split("Format:")[1].strip().replace("\"","").split("|")

    ## Get the ConSplice info and create an index
    consplice_info_index = vcf["ConSplice"]["Description"].strip().split("Format:")[1].strip().replace("\"","").split("|")

    #### The SQUIRLS annotation does not have Format info

    ## Add ConSplice to header
    vcf.add_info_to_header({"ID": "ConSpliceML", "Description": "The ConSpliceML pathogenic splice-altering prediction. Format: Gene_Name|ConSpliceML_prediction", "Type": "String", "Number": "."})


    ## Open vcf writer
    scored_vcf = (Writer(outfile, vcf, "w") if args.out_type == "vcf" 
                  else Writer(outfile, vcf, "wz") if args.out_type == "vcfgz"
                  else Writer(outfile, vcf, "wbu") if args.out_type == "bcf"
                  else Writer(outfile, vcf, "wb") if args.out_type == "bcfgz"
                  else Writer(outfile, vcf, "w"))


    ## Iterate over variants
    print("\nParsing vcf file and using ConSpliceML to score variants")
    for var in vcf:
        
        ##Don't try to score variants that are missing a SpliceAI score, a SQUIRLS score, or a ConSplice score
        if var.INFO["SpliceAI"] is None or var.INFO["SQUIRLS_SCORE"] is None or var.INFO["ConSplice"] is None:
            scored_vcf.write_record(var)
            continue


        ## Get the SpliceAI annotation as a dict of dicts
        spliceai_annotation = {dict(zip(spliceai_info_index, x.strip().split("|")))["SYMBOL"]:dict(zip(spliceai_info_index, x.strip().split("|"))) for x in var.INFO["SpliceAI"].strip().split(",")}


        ## Get MAX SpliceAI score per annotation
        for gene, anno in spliceai_annotation.items():
            
            spliceai_annotation[gene]["MAX_SCORE"] = max([float(anno["DS_AG"]), 
                                                          float(anno["DS_AL"]), 
                                                          float(anno["DS_DG"]), 
                                                          float(anno["DS_DL"])])


        ## Get the ConSplice annotation as a dict of dicts
        consplice_annotation = {dict(zip(consplice_info_index, x.strip().split("|")))["Gene_Name"]:dict(zip(consplice_info_index, x.strip().split("|"))) for x in var.INFO["ConSplice"].strip().split(",")}


        ## Get the Max SQUIRLS score (SQUIRLS does not have a defined format other than a score for each transcript.)
        max_squirls = max([float(x.strip().split("=")[1]) for x in var.INFO["SQUIRLS_SCORE"].strip().split("|")[1:]])


        ## Match SpliceAI and ConSplice scores by gene name 
        combined_score_list = []

        ### Iterate over each ConSplice gene annotation
        for gene_name in consplice_annotation:
            
            final_gene_name = None  

            ### Check for matching gene symbol
            if gene_name in spliceai_annotation:

                final_gene_name = gene_name

            ### Check for matching alt gene symobl
            else:
                
                for gene in alt_symbol_dict[gene_name]:

                    if gene in spliceai_annotation:
                        
                        final_gene_name = gene

                        break
                    
            if final_gene_name is None:
                continue
                
            ## Add scores
            combined_score_list.append({"Gene_Name": gene_name,
                                        metadata["Features"]["SpliceAI"]:  spliceai_annotation[final_gene_name]["MAX_SCORE"], 
                                        metadata["Features"]["SQUIRLS"]:   max_squirls, 
                                        metadata["Features"]["ConSplice"]: float(consplice_annotation[gene_name]["ConSplice_Score"])}) 


        ## Don't score this variant if there is no matching gene symbols between SpliceAI and ConSplice
        if len(combined_score_list) < 1:
            scored_vcf.write_record(var)
            continue

        ## Create a variant dataframe with feature score columns
        var_df = pd.DataFrame(combined_score_list)
        var_df = var_df[(["Gene_Name"] + feature_order)]

        ## Score using the ConSpliceML trained model
        scored_var_df, ConSpliceML_col  = ConSpliceML_score(rf = rf,
                                                            var_df = var_df,
                                                            feature_col_names = feature_order)

        ## Create the new annotation
        new_anno = ["{}|{}".format(d["Gene_Name"],d[ConSpliceML_col]) for d in scored_var_df[["Gene_Name",ConSpliceML_col]].to_dict("records")]

        ## Add the new annotation to the variant  
        var.INFO["ConSpliceML"] = ",".join(new_anno)
        
        ## Write variant to output vcf file
        scored_vcf.write_record(var)


    ## Close cyvcf2 IO objects
    vcf.close()
    scored_vcf.close()

    print("\nVariants in the VCF file have been scored by ConSpliceML and saved as '{}'".format(outfile))
        
    print("\nDONE")

