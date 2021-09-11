import argparse
import sys

from .__init__ import __version__
from .ML.rf_training import add_training
from .ML.rf_scoring import add_ml_scoring
from .constarint.create_substitution_matrix import add_substitution_matrix
from .constarint.regional_o_e_scores import add_regional_o_e
from .constarint.gene_o_e_scores import add_gene_o_e
from .constarint.calculate_constraint_score import add_constraint_scores
from .constarint.agg_regional_scores import add_aggregate_overlapping_regions
from .constarint.select_scores import add_select_score
from .constarint.convert_to_bed import add_to_bed



if sys.version_info[0] < 3:
    print("[ConSplice_cli:ERROR] Python 2 is not supported. ConSplice requires python 3")
    sys.exit(1)


def main(args=None):
    if args is None:
            args = sys.argv[1:]
    
    parser = argparse.ArgumentParser(
        prog="ConSplice", 
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description = ("\n\t===================================================\n"
                       "\t| ConSplice-CLI: The CLI to the ConSplice project |\n"
                       "\t===================================================\n\n "
                       "\t(1) Score/Annotate variants using the ConSpliceML model to identify potential pathogenic clinically relevant alternative splicing variants." 
                       "\n\t(2) Train a Random Forest model for deleterious splicing prediction and interpritation utilizing the constrained splicing model."
                       "\n\t(3) Generate genic or regional splicing constraint profiles.") 


    )


    parser.add_argument(
        "-v",
        "--version",
        help = "Installed ConSplice CLI Version",
        action="version",
         version="%(prog)s v" + str(__version__),
    )

    ## 1st lay sub command
    sub = parser.add_subparsers(title="Main Commands", dest="command")
    sub.required = True

    
    ## 2nd layer sub command for ConSplice ML function
    ml = sub.add_parser("ML", 
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                        help = "ConSpliceML - Score variants with ConSpliceML or train a new ConSpliceML model", 
                        description = ("\n\t###############\n"
                                       "\t# ConSpliceML #\n"
                                       "\t###############\n\n"
                                       "\tModule to:\n"
                                       "\t - Score variants to deterime pathogenic alterniative splicing variants using the ConSpliceML model.\n"
                                       "\t - Train a new ConSpliceML model."
                        )
    )

    ## 2nd layer sub command for Constraint function
    con = sub.add_parser("constraint", 
                        formatter_class=argparse.RawDescriptionHelpFormatter,
                         help = "Generate genic or regional splicing constraint profiles", 
                         description = ("\n\t#######################\n"
                                        "\t# Splicing Constraint #\n"
                                        "\t#######################\n\n"
                                        "\tModule to generate genic or regional\n"
                                        "\tsplicing constraint using patterns of\n"
                                        "\tpurifying selection and evidence of\n"
                                        "\talternative splicing"
                        )
    )



    ## ML functions
    sub_ml = ml.add_subparsers()

    add_training(sub_ml)

    add_ml_scoring(sub_ml)



    ## Constraint functions
    sub_con = con.add_subparsers()

    add_substitution_matrix(sub_con)

    add_regional_o_e(sub_con)

    add_gene_o_e(sub_con)

    add_constraint_scores(sub_con)

    add_aggregate_overlapping_regions(sub_con)

    add_select_score(sub_con)

    add_to_bed(sub_con)



    ## Initate argparse
    args = parser.parse_args(args)
    args.func(parser, args)




if __name__ == "__main__":
    sys.exit(main() or 0)
