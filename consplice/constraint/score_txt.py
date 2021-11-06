from __future__ import print_function
import sys
import io
import os
import argparse 
from .utils import extract_constraint_score, get_alternative_gene_symbols


#---------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
#---------------------------------------------------------------------------------------------------------------------------------

def add_score_txt(sub_p):

    p = sub_p.add_parser("score-txt",
                         help = "Add ConSplice scores to a tab-delimited txt variant file",
                         formatter_class=argparse.RawDescriptionHelpFormatter,
                         description=("\n\t********************************************\n"
                                      "\t* ConSplice - Score tab-delimited txt with *\n"
                                      "\t********************************************\n\n"
                                      "\tScore a 1-based variant file in tab-delimited txt format with ConSplice scores\n"
                                      "\t - Expecting a tab delimited file with a header line that start with a #\n"
                                      "\t - The last header line with a # will be treated as a column index based\n"
                                      "\t    on the tab seperated names in that line.\n"
                                      "\t - The variant position will be treated as 1-based genomic position\n"
                                      "\t - The ConSplice score will be assigned by either a gene name or the \n"
                                      "\t    max ConSplice score for that position based on input parameters.\n"

                         )
    )

    req = p.add_argument_group("Required Arguments")
    rif = p.add_argument_group("Required if using 'by-gene'")
    req.add_argument(
        "--consplice-file",
        metavar = "ConSplice Score bed File", 
        required=True, 
        help="(Required) The path to the 0-based ConSplice score file in bed format."
    )

    req.add_argument(
        "--txt-file",
        metavar = "Variant txt file", 
        required=True, 
        help="(Required) The path to the 1-based variant file in txt format."
    )

    req.add_argument(
        "--score-type",
        choices = ["by-gene","max"],
        help = "(Required) How to add the score to the variant file. Either by the gene name for the variant or the max ConSplice score for that position. Choices = 'by-gene' or 'max') (NOTE: If 'by-gene is chosen, then the --txt-gene-name argument will be used to determine the gene name in the txt file)"
    )

    req.add_argument(
        "--out-file",
        metavar="Ouput File",
        required = True,
        help = "(Required) The path and/or the name of the output file to create"
    )

    p.add_argument(
        "--txt-chrom",
        metavar = "txt chromosome column", 
        default = "chrom",
        help="The name of the chromosome column in the txt file. (Default = chrom) (NOTE: Expecting tab delimited file with a header line)"
    )

    p.add_argument(
        "--txt-pos",
        metavar = "txt position column", 
        default = "pos",
        help="The name of the position column in the txt file. (Default = pos) (NOTE: Expecting tab delimited file with a header line)"
    )

    p.add_argument(
        "--txt-gene-name",
        metavar = "txt gene name column", 
        default = "gene_name",
        help="The name for the gene name/symbol column in the txt file. (Default = gene_name) (NOTE: Expecting tab delimited file with a header line)"
    )

    rif.add_argument(
        "--alt-gene-symbol",
        metavar = "Alternative Gene Symbol File",
        default = "NA",
        help = "(Required if --score-type == 'by-gene') A file that contains mappings between a canonical gene symbol to alternative gene symbols. NOTE: The file needs to have a header!. This script is set up to use the HGNC protein-coding gene mapping file can be found in the ConSplice GitHub repo: https://github.com/mikecormier/ConSplice"  
    )

    p.add_argument(
        "--consplice-col",
        metavar = "ConSplice Column Name",
        default = "ConSplice_Percentile",
        help = "The name of the ConSplice score column in the ConSplice file. (Default = ConSplice_Percentile)"
    )

    p.add_argument(
        "--out-consplice-col",
        metavar="ConSplice Score Column",
        default = "ConSplice_score",
        help = "The name of the column to create the represents the ConSplice score. (Default = 'ConSplice_score')"
    )

    p.set_defaults(func=add_conSplice_score)






#---------------------------------------------------------------------------------------------------------------------------------
## Main
#---------------------------------------------------------------------------------------------------------------------------------


def add_conSplice_score(parser, args):

    print("\n\n\t***************************************")
    print("\t* ConSplice - Score tab-delimited txt *")
    print("\t***************************************\n\n")


    print(("\nInput Arguments:"
           "\n================"
           "\n - consplice-file:           {}"
           "\n - txt-file:                 {}"
           "\n - score-type:               {}"
           "\n - out-file:                 {}"
           "\n - txt-chrom:                {}" 
           "\n - txt-pos:                  {}"
           "\n - txt-gene-name:            {}"
           "\n - alt_gene_symbol:          {}"
           "\n - consplice-col:            {}"
           "\n - out-consplice-col:        {}"
           "\n"
           ).format(args.consplice_file,
                    args.txt_file, 
                    args.score_type,
                    args.out_file,
                    args.txt_chrom,
                    args.txt_pos,
                    args.txt_gene_name,
                    args.alt_gene_symbol,
                    args.consplice_col,
                    args.out_consplice_col,
                    )
    )


    ## Check file paths
    print("\nChecking file paths")
    assert os.path.exists(args.consplice_file), "\n!!ERROR!! The ConSplice Score file does not exists"
    assert os.path.exists(args.txt_file), "\n!!ERROR!! The variant txt file does not exists"

    if args.score_type == "by-gene":
        assert os.path.exists(args.alt_gene_symbol), "\n!!ERROR!! The alternative gene symbols file does not exists"
        
        print("\nCreating an map of alternative gene symbols")
        ## Get a dictionary of alternative gene symbols
        alt_symbol_dict = get_alternative_gene_symbols(args.alt_gene_symbol)


    print("\nChecking for header in variant file")
    header = []
    h_index = 0
    with io.open(args.txt_file, "rt", encoding = "utf-8") as fh:
        
        header_line = fh.readline()

        assert header_line[0] == "#", "\n!!ERROR!! the first line of the txt file should be a header line with a # starting the header" 

        header = header_line.replace("#","").strip().split("\t")

        ## Iterate over the file until the header lines starting with a # are done
        for line in fh:
            
            if line[0] == "#":
                header = line.replace("#","").strip().split("\t")
            else:
                break

            ## update header index
            h_index += 1

        assert args.txt_chrom in header, "\n!!ERROR!! the '{}' chrom column is not in the txt file".format(args.txt_chrom)
        assert args.txt_pos in header, "\n!!ERROR!! the '{}' position column is not in the txt file".format(args.txt_pos)

        if args.score_type == "by-gene":
            assert args.txt_gene_name in header, "\n!!ERROR!! the '{}' gene name column is not in the txt file".format(args.txt_gene_name)


    print("\nLoading ConSplice scores into an interval tree")

    consplice_interlap = extract_constraint_score(file_path  = args.consplice_file, 
                                                  chrom_col  = "chrom",
                                                  start_col  = "region_start",
                                                  end_col    = "region_end",
                                                  score_col  = args.consplice_col,
                                                  zero_based = True if args.consplice_file.endswith(".bed") or args.consplice_file.endswith("bed.gz") else False)



    print("\nParsing variant file and adding ConSplice scores")
    with io.open(args.txt_file, "rt", encoding = "utf-8") as fh, io.open(args.out_file, "w") as out_fh:
        
        ## iterate over each variant in the file
        for i,line in enumerate(fh):
            
            ## Write the header lines 
            if line[0] == "#":
                line = line.strip() + "\t{}\n".format(args.out_consplice_col) if i == h_index else line
                out_fh.write(line)
                continue
            
            line_list = line.strip().split("\t")
            line_dict = dict(zip(header, line_list))

            chrom = line_dict[args.txt_chrom].replace("chr","")
            pos = int(line_dict[args.txt_pos])

            ## Get ConSplice score
            var_score = float("NaN")

            consplice_scores = list(consplice_interlap[str(chrom)].find((pos,pos)))

            if args.score_type == "max":
                var_score = max([float(x[2]) for x in consplice_scores])

            elif args.score_type == "by-gene":
                
                gene_name = line_dict[args.txt_gene_name].strip()

                score_list = []
                for score in consplice_scores:
                    
                    ## Update the alt_symbol_dict with the current gene. (Some genes may not be included in the alt symbol file) 
                    alt_symbol_dict[score[3]].add(score[3])

                    if gene_name in alt_symbol_dict[score[3]]:
                        
                        score_list.append(score[2])


                var_score = max(score_list) if len(score_list) > 0 else var_score 

            line_list.append(str(var_score))

            out_fh.write("\t".join(line_list) + "\n")


    print("\nConSplice scores added to '{}' output file under the '{}' column name".format(args.out_file, args.out_consplice_col))

    print("\nDone")
