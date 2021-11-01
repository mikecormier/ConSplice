from __future__ import print_function
import sys
import io
import os
import argparse 
from collections import defaultdict
from cyvcf2 import VCF,Writer
from .utils import extract_constraint_score


#---------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
#---------------------------------------------------------------------------------------------------------------------------------

def add_score_vcf(sub_p):

    p = sub_p.add_parser("score-vcf",
                         help = "Add ConSplice scores to a vcf file",
                         formatter_class=argparse.RawDescriptionHelpFormatter,
                         description=("\n\t****************************\n"
                                      "\t* Score vcf with ConSplice *\n"
                                      "\t****************************\n\n"
                                      "\tScore a valid vcf file with ConSplice scores\n"
                                      "\t - Each variant is treated as an SNV using the\n"
                                      "\t   variant position rather the the genomic range\n"

                         )
    )

    req = p.add_argument_group("Required Arguments")
    req.add_argument(
        "--consplice-file",
        metavar = "ConSplice Score bed File", 
        required=True, 
        help="(Required) The path to the 0-based ConSplice score file in bed format."
    )

    req.add_argument(
        "--vcf-file",
        metavar = "The VCF file", 
        required=True, 
        help="(Required) The path to the valid vcf file to add ConSplice scores to"
    )

    req.add_argument(
        "--out-file",
        metavar="Ouput File",
        required = True,
        help = "(Required) The path and/or the name of the output file to create"
    )
    
    p.add_argument(
        "--out-type",
        metavar="Output file type",
        choices = ["vcf","vcfgz","bcf","bcfgz"],
        default = "vcf",
        help="The output file type. Whether the output file should be a normal vcf/bcf file or a bgzipped vcf/bcf file. Choices = ['vcf','vcfgz','bcf','bcfgz']. vcfgz and bcfgz are the compressed version of vcf or bcf files. Default = 'vcf'"
    )

    p.add_argument(
        "--consplice-col",
        metavar = "ConSplice Column Name",
        default = "ConSplice_percentile",
        help = "The name of the ConSplice score column in the ConSplice file. (Default = ConSplice_percentile)"
    )

    p.set_defaults(func=add_conSplice_score)






#---------------------------------------------------------------------------------------------------------------------------------
## Main
#---------------------------------------------------------------------------------------------------------------------------------


def add_conSplice_score(parser, args):

    print("\n\n\t****************************")
    print("\t* Score vcf with ConSplice *")
    print("\t****************************\n\n")



    outfile = args.out_file
    if not args.out_file.endswith("vcf") and args.out_type == "vcf":
        outfile = args.out_file + ".vcf"

    elif not args.out_file.endswith("bcf") and args.out_type == "bcf":
        outfile = args.out_file + ".bcf"

    elif not args.out_file.endswith("vcf.gz") and args.out_type == "vcfgz":
        outfile = args.out_file + ".vcf.gz" if not args.out_file.endswith(".vcf") else args.out_file + ".gz"

    elif not args.out_file.endswith("bcf") and args.out_type == "bcfgz":
        outfile = args.out_file + ".bcf" 

    print(("\nInput Arguments:"
           "\n================"
           "\n - consplice-file:           {}"
           "\n - vcf-file:                 {}"
           "\n - out-file:                 {}"
           "\n - out-type:                 {}" 
           "\n - consplice-col:            {}"
           "\n"
           ).format(args.consplice_file,
                    args.vcf_file, 
                    outfile,
                    args.out_type,
                    args.consplice_col,
                    )
    )


    ## Check file paths
    print("\nChecking file paths")
    assert os.path.exists(args.consplice_file), "\n!!ERROR!! The ConSplice Score file does not exists"
    assert os.path.exists(args.vcf_file), "\n!!ERROR!! The vcf file does not exists"


    print("\nLoading ConSplice scores into an interval tree")

    consplice_interlap = extract_constraint_score(file_path  = args.consplice_file, 
                                                  chrom_col  = "chrom",
                                                  start_col  = "region_start",
                                                  end_col    = "region_end",
                                                  score_col  = args.consplice_col,
                                                  zero_based = True if args.consplice_file.endswith(".bed") or args.consplice_file.endswith("bed.gz") else False)


    print("\nParsing vcf file and adding ConSplice scores")

    ## Open vcf file
    vcf = VCF(args.vcf_file)

    ## Add ConSplice to header
    vcf.add_info_to_header({"ID": "ConSplice", "Description": "The ConSplice percentile score annotation. Format: Gene_Name|ConSplice_Score", "Type": "Float", "Number": "."})

    ## Open vcf writer
    scored_vcf = (Writer(outfile, vcf, "w") if args.out_type == "vcf" 
                  else Writer(outfile, vcf, "wz") if args.out_type == "vcfgz"
                  else Writer(outfile, vcf, "wbu") if args.out_type == "bcf"
                  else Writer(outfile, vcf, "wb") if args.out_type == "bcfgz"
                  else Writer(outfile, vcf, "w"))


    ## Iterate over variants
    for var in vcf:
        
        ## Get chrom and position of the variant
        chrom = str(var.CHROM).replace("chr","")
        pos = int(var.POS)

        ## Get the ConSplice scores for the current variant 
        score_dict = defaultdict(list)
        for score in consplice_interlap[str(chrom)].find((pos,pos)):
            
            score_dict[score[3]].append(score[2])

        score_list = ["{}|{}".format(gene, max(scores)) for gene,scores in score_dict.items()]

        if score_list:
            var.INFO["ConSplice"] = ",".join(score_list)

        ## Write variant to output vcf file
        scored_vcf.write_record(var)

    ## Close cyvcf2 IO objects
    vcf.close()
    scored_vcf.close()


    print("\nConSplice scores added to '{}' output file using the annotation name 'ConSplice'".format(outfile))

    print("\nDone")
