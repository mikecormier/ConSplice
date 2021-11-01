import argparse
import sys
import os

from .__init__ import __version__, __cur_path__
from .ML.rf_training import add_training
from .ML.rf_scoring import add_ml_scoring
from .constraint.create_substitution_matrix import add_substitution_matrix
from .constraint.o_e_counts import add_o_e_counts
from .constraint.calculate_constraint_score import add_constraint_scores
from .constraint.agg_regional_scores import add_aggregate_overlapping_regions
from .constraint.select_scores import add_select_score
from .constraint.convert_to_bed import add_to_bed
from .constraint.score_txt import add_score_txt
from .constraint.score_bed import add_score_bed
from .constraint.score_vcf import add_score_vcf




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

    parser.add_argument(
        "--config-path",
        default = os.path.join(os.path.dirname(str(__cur_path__)), "../config/config.yml"),
        help = "(constraint subcommand only) File path to the new ConSplice config yaml file. By default, ConSplice CLI will use the default ConSplice config.yml file",
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

    add_o_e_counts(sub_con)

    add_constraint_scores(sub_con)

    add_aggregate_overlapping_regions(sub_con)

    add_select_score(sub_con)

    add_to_bed(sub_con)

    add_score_txt(sub_con)

    add_score_bed(sub_con)

    add_score_vcf(sub_con)



    ## Initate argparse
    args = parser.parse_args(args)
    args.func(parser, args)




if __name__ == "__main__":
    sys.exit(main() or 0)
