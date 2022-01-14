import argparse
import io
import os
import sys

from Bio.bgzf import BgzfWriter

# ---------------------------------------------------------------------------------------------------------------------------------
## Argument Parser
# ---------------------------------------------------------------------------------------------------------------------------------


def add_to_bed(sub_p):

    p = sub_p.add_parser(
        "to-bed",
        help="Convert the 1-based scored ConSplice txt file to a 0-based bed file",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "\n\t**********************\n"
            "\t* ConSplice - to bed *\n"
            "\t**********************\n\n"
            "\tConvert the scored ConSplice txt file, where the genomic positions are\n"
            "\t1-based like a vcf file, to a bed file, where the genomic positions will\n"
            "\tbe 0-based"
        ),
    )

    req = p.add_argument_group("Required Arguments")
    req.add_argument(
        "--score-file",
        metavar="Regional Score File",
        required=True,
        help="(Required) The path to the scored ConSplice txt file to convert to 0-based bed file",
    )

    req.add_argument(
        "--out-file",
        metavar="Output file",
        required=True,
        help="(Required) The name of the output bed file to create.",
    )

    p.add_argument(
        "--out-type",
        metavar="Output file type",
        choices=["bed", "bedgz"],
        default="bed",
        help="(Optional) The output file type. Whether the output file should be a normal bed file or a bgzipped bed file. Choices = 'bed' or 'bedgz'. Default = 'bed'",
    )

    p.set_defaults(func=to_bed)


# ---------------------------------------------------------------------------------------------------------------------------------
## Main
# ---------------------------------------------------------------------------------------------------------------------------------


def to_bed(parser, args):

    outfile = args.out_file
    if not args.out_file.endswith("bed") and args.out_type == "bed":
        outfile = args.out_file + ".bed"

    elif not args.out_file.endswith("bed.gz") and args.out_type == "bedgz":
        outfile = (
            args.out_file + ".bed.gz"
            if not args.out_file.endswith(".bed")
            else args.out_file + ".gz"
        )

    print("\n\t**********************")
    print("\t* ConSplice - to bed *")
    print("\t**********************\n\n")

    print(
        (
            "\nInput Arguments:"
            "\n================"
            "\n - score-file:           {}"
            "\n - out-file:             {}"
            "\n - out-type:             {}"
        ).format(
            args.score_file,
            outfile,
            args.out_type,
        )
    )

    print(
        "\nParsing input file, updating to 0-based and writing to '{}'".format(outfile)
    )

    ## open ConSplice file
    with io.open(args.score_file, "rt", encoding="utf-8") as in_fh:

        ## Open output file
        try:
            out_fh = (
                open(outfile, "w")
                if args.out_type == "bed"
                else BgzfWriter(outfile, "wb")
            )

        except IOError as e:
            print(str(e))
            sys.exit(1)

        ## Update output file with header
        header_line = in_fh.readline()

        out_fh.write(header_line)

        header = header_line.strip().split("\t")

        header_dict = dict(zip(header, range(len(header))))

        ## Add each line to the output file with the start position adjusted from 1-based to 0-based
        for line in in_fh:

            line_list = line.strip().split("\t")

            line_list[header_dict["region_start"]] = str(
                int(line_list[header_dict["region_start"]]) - 1
            )

            out_fh.write("\t".join(line_list) + "\n")

        out_fh.close()

    print("\nNew file written to '{}'".format(outfile))

    print("\nDONE")
