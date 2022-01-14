from __future__ import print_function

import argparse
import gzip
import io
import os
import sys

from Bio.bgzf import BgzfWriter

from .utils import extract_constraint_score, get_alternative_gene_symbols

# ---------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
# ---------------------------------------------------------------------------------------------------------------------------------


def add_score_bed(sub_p):

    p = sub_p.add_parser(
        "score-bed",
        help="Add ConSplice scores to a bed variant file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "\n\t*************************\n"
            "\t* ConSplice - Score bed *\n"
            "\t*************************\n\n"
            "\tScore a 0-based variant file in bed format with ConSplice scores\n"
            "\t - Expecting a tab delimited file with a header line that start with a #\n"
            "\t - The last header line with a # will be treated as a column index based\n"
            "\t    on the tab separated names in that line.\n"
            "\t - The variant position will be treated as 0-based genomic position\n"
            "\t - The ConSplice score will be assigned by either a gene name or the \n"
            "\t    max ConSplice score for that position based on input parameters.\n"
        ),
    )

    req = p.add_argument_group("Required Arguments")
    rif = p.add_argument_group("Required if using 'by-gene'")
    req.add_argument(
        "--consplice-file",
        metavar="ConSplice Score bed File",
        required=True,
        help="(Required) The path to the 0-based ConSplice score file in bed format.",
    )

    req.add_argument(
        "--bed-file",
        metavar="Variant bed file",
        required=True,
        help="(Required) The path to the 0-based variant file in bed format.",
    )

    req.add_argument(
        "--score-type",
        choices=["by-gene", "max"],
        help="(Required) How to add the score to the variant file. Either by the gene name for the variant or the max ConSplice score for that position. Choices = 'by-gene' or 'max') (NOTE: If 'by-gene is chosen, then the --bed-gene-name argument will be used to determine the gene name in the bed file)",
    )

    req.add_argument(
        "--out-file",
        metavar="Ouput File",
        required=True,
        help="(Required) The path and/or the name of the output file to create",
    )

    p.add_argument(
        "--out-type",
        metavar="Output file type",
        choices=["bed", "bedgz"],
        default="bed",
        help="The output file type. Whether the output file should be a normal bed file or a bgzipped bed file. Choices = ['bed', 'bedgz']. Default = 'bed'",
    )

    p.add_argument(
        "--bed-chrom",
        metavar="bed chromosome column",
        default="chrom",
        help="The name of the chromosome column in the bed file. (Default = chrom) (NOTE: Expecting tab delimited file with a header line)",
    )

    p.add_argument(
        "--bed-pos",
        metavar="bed start position column",
        default="start",
        help="The name of the start position column in the bed file. (Default = start) (NOTE: Expecting tab delimited file with a header line)",
    )

    p.add_argument(
        "--bed-gene-name",
        metavar="bed gene name column",
        default="gene_name",
        help="The name for the gene name/symbol column in the bed file. (Default = gene_name) (NOTE: Expecting tab delimited file with a header line)",
    )

    rif.add_argument(
        "--alt-gene-symbol",
        metavar="Alternative Gene Symbol File",
        default="NA",
        help="(Required if --score-type == 'by-gene') A file that contains mappings between a canonical gene symbol to alternative gene symbols. NOTE: The file needs to have a header!. This script is set up to use the HGNC protein-coding gene mapping file can be found in the ConSplice GitHub repo: https://github.com/mikecormier/ConSplice",
    )

    p.add_argument(
        "--consplice-col",
        metavar="ConSplice Column Name",
        default="ConSplice_Percentile",
        help="The name of the ConSplice score column in the ConSplice file. (Default = ConSplice_Percentile)",
    )

    p.add_argument(
        "--out-consplice-col",
        metavar="ConSplice Score Column",
        default="ConSplice_score",
        help="The name of the column to create the represents the ConSplice score. (Default = 'ConSplice_score')",
    )

    p.set_defaults(func=add_conSplice_score)


# ---------------------------------------------------------------------------------------------------------------------------------
## Main
# ---------------------------------------------------------------------------------------------------------------------------------


def add_conSplice_score(parser, args):

    print("\n\n\t*************************")
    print("\t* ConSplice - Score bed *")
    print("\t*************************\n\n")

    outfile = args.out_file
    if not args.out_file.endswith("bed") and args.out_type == "bed":
        outfile = args.out_file + ".bed"

    elif not args.out_file.endswith("bed.gz") and args.out_type == "bedgz":
        outfile = (
            args.out_file + ".bed.gz"
            if not args.out_file.endswith(".bed")
            else args.out_file + ".gz"
        )

    print(
        (
            "\nInput Arguments:"
            "\n================"
            "\n - consplice-file:           {}"
            "\n - bed-file:                 {}"
            "\n - score-type:               {}"
            "\n - out-file:                 {}"
            "\n - out-type:                 {}"
            "\n - bed-chrom:                {}"
            "\n - bed-pos:                  {}"
            "\n - bed-gene-name:            {}"
            "\n - alt_gene_symbol:          {}"
            "\n - consplice-col:            {}"
            "\n - out-consplice-col:        {}"
            "\n"
        ).format(
            args.consplice_file,
            args.bed_file,
            args.score_type,
            outfile,
            args.out_type,
            args.bed_chrom,
            args.bed_pos,
            args.bed_gene_name,
            args.alt_gene_symbol,
            args.consplice_col,
            args.out_consplice_col,
        )
    )

    ## Check file paths
    print("\nChecking file paths")
    assert os.path.exists(
        args.consplice_file
    ), "\n!!ERROR!! The ConSplice Score file does not exists"
    assert os.path.exists(
        args.bed_file
    ), "\n!!ERROR!! The variant bed file does not exists"

    assert args.consplice_file.endswith(".bed") or args.consplice_file.endswith(
        ".bed.gz"
    ), "\n!!ERROR!! Expecting a .bed or .bed.gz ConSplice file"

    if args.score_type == "by-gene":
        assert os.path.exists(
            args.alt_gene_symbol
        ), "\n!!ERROR!! The alternative gene symbols file does not exists"

        print("\nCreating an map of alternative gene symbols")
        ## Get a dictionary of alternative gene symbols
        alt_symbol_dict = get_alternative_gene_symbols(args.alt_gene_symbol)

    if os.path.exists(outfile):
        print(
            "\n\n**WARNING** The output file '{}' already exists. This file will be overwritten. To halt this process press CTRL+C\n".format(
                outfile
            )
        )

    print("\nChecking for header in variant file")
    header = []
    h_index = 0

    try:
        fh = (
            gzip.open(args.bed_file, "rt", encoding="utf-8")
            if args.bed_file.endswith(".gz")
            else io.open(args.bed_file, "rt", encoding="utf-8")
        )
    except IOError as e:
        print("\n!!ERROR!! Unable to read bed file\n")
        print(str(e))
        sys.exit(1)

    header_line = fh.readline()

    assert (
        header_line[0] == "#"
    ), "\n!!ERROR!! the first line of the bed file should be a header line with a # starting the header"

    header = header_line.replace("#", "").strip().split("\t")

    ## Iterate over the file until the header lines starting with a # are done
    for line in fh:

        if line[0] == "#":
            header = line.replace("#", "").strip().split("\t")
        else:
            break

        ## update header index
        h_index += 1

    fh.close()

    assert (
        args.bed_chrom in header
    ), "\n!!ERROR!! the '{}' chrom column is not in the bed file".format(args.bed_chrom)
    assert (
        args.bed_pos in header
    ), "\n!!ERROR!! the '{}' position column is not in the bed file".format(
        args.bed_pos
    )

    if args.score_type == "by-gene":
        assert (
            args.bed_gene_name in header
        ), "\n!!ERROR!! the '{}' gene name column is not in the bed file".format(
            args.bed_gene_name
        )

    print("\nLoading ConSplice scores into an interval tree")

    consplice_interlap = extract_constraint_score(
        file_path=args.consplice_file,
        chrom_col="chrom",
        start_col="region_start",
        end_col="region_end",
        score_col=args.consplice_col,
        zero_based=False,
    )  ## Don't adjust for for zero based

    print("\nParsing variant file and adding ConSplice scores")

    ## Open files
    try:
        fh = (
            gzip.open(args.bed_file, "rt", encoding="utf-8")
            if args.bed_file.endswith(".gz")
            else io.open(args.bed_file, "rt", encoding="utf-8")
        )
        out_fh = (
            io.open(outfile, "w")
            if args.out_type == "bed"
            else BgzfWriter(outfile, "wb")
        )
    except IOError as e:
        print("\n!!ERROR!! Unable to open file\n")
        print(str(e))
        sys.exit(1)

    ## iterate over each variant in the file
    for i, line in enumerate(fh):

        ## Write the header lines
        if line[0] == "#":
            line = (
                line.strip() + "\t{}\n".format(args.out_consplice_col)
                if i == h_index
                else line
            )
            out_fh.write(line)
            continue

        line_list = line.strip().split("\t")
        line_dict = dict(zip(header, line_list))

        chrom = line_dict[args.bed_chrom].replace("chr", "")
        pos = int(line_dict[args.bed_pos])

        ## Get ConSplice score
        var_score = float("NaN")

        consplice_scores = list(consplice_interlap[str(chrom)].find((pos, pos)))

        if args.score_type == "max":
            var_score = max([float(x[2]) for x in consplice_scores])

        elif args.score_type == "by-gene":

            gene_name = line_dict[args.bed_gene_name].strip()

            score_list = []
            for score in consplice_scores:

                if gene_name in alt_symbol_dict[score[3]]:

                    score_list.append(score[2])

            var_score = max(score_list) if len(score_list) > 0 else var_score

        line_list.append(str(var_score))

        out_fh.write("\t".join(line_list) + "\n")

    fh.close()
    out_fh.close()

    print(
        "\nConSplice scores added to '{}' output file under the '{}' column name".format(
            outfile, args.out_consplice_col
        )
    )

    print("\nDone")
