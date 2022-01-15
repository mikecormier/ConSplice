import argparse
import os
import sys

from .__init__ import __cur_path__, __version__
from .constraint.agg_regional_scores import add_aggregate_overlapping_regions
from .constraint.calculate_constraint_score import add_constraint_scores
from .constraint.convert_to_bed import add_to_bed
from .constraint.create_substitution_matrix import add_substitution_matrix
from .constraint.o_e_counts import add_o_e_counts
from .constraint.score_bed import add_score_bed
from .constraint.score_txt import add_score_txt
from .constraint.score_vcf import add_score_vcf
from .constraint.select_scores import add_select_score
from .ML.rf_scoring_vcf import add_ml_scoring_vcf
from .ML.rf_training import add_training

if sys.version_info[0] < 3:
    print(
        "[ConSplice_cli:ERROR] Python 2 is not supported. ConSplice requires python 3"
    )
    sys.exit(1)


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="ConSplice",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "\n\t=====================================================\n"
            "\t|| ConSplice-CLI: The CLI to the ConSplice project ||\n"
            "\t=====================================================\n\n"
            "\t(1) Score/Annotate variants using the ConSpliceML model to identify potential pathogenic clinically relevant alternative splicing variants."
            "\n\t(2) Train a Random Forest model for deleterious splicing prediction and interpretation utilizing the constrained splicing model."
            "\n\t(3) Generate genic or regional splicing constraint profiles."
        ),
    )

    parser.add_argument(
        "-v",
        "--version",
        help="Installed ConSplice CLI Version",
        action="version",
        version="%(prog)s v" + str(__version__),
    )

    parser.add_argument(
        "--config-path",
        default=os.path.join(
            os.path.dirname(str(__cur_path__)), "config/config.yml"
        ),
        help="(constraint subcommand only) File path to the new ConSplice config yaml file. By default, ConSplice CLI will use the default ConSplice config.yml file",
    )

    parser.add_argument(
        "--base-config",
        default=os.path.join(os.path.dirname(str(__cur_path__)), "config/"),
        help=argparse.SUPPRESS,
    )

    parser.add_argument(
        "--check-config-path",
        action="store_true",
        help=argparse.SUPPRESS
    )


    ## 1st lay sub command
    sub = parser.add_subparsers(title="Main Commands", dest="command")
    #sub.required = True

    ## 2nd layer sub command for ConSplice ML function
    ml = sub.add_parser(
        "ML",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="ConSpliceML - Score variants with ConSpliceML or train a new ConSpliceML model",
        description=(
            "\n\t###############\n"
            "\t# ConSpliceML #\n"
            "\t###############\n\n"
            "\tModule to:\n"
            "\t - Score variants to determine pathogenic alternative splicing variants using the ConSpliceML model.\n"
            "\t - Train a new ConSpliceML model."
        ),
    )

    ## 2nd layer sub command for Constraint function
    con = sub.add_parser(
        "constraint",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Generate genic or regional splicing constraint profiles",
        description=(
            "\n\t#######################\n"
            "\t# Splicing Constraint #\n"
            "\t#######################\n\n"
            "\tModule to generate genic or regional\n"
            "\tsplicing constraint using patterns of\n"
            "\tpurifying selection and evidence of\n"
            "\talternative splicing"
        ),
    )

    ## ML functions
    sub_ml = ml.add_subparsers()

    add_training(sub_ml)

    add_ml_scoring_vcf(sub_ml)

    ## Constraint functions
    sub_con = con.add_subparsers()

    add_substitution_matrix(sub_con)

    add_o_e_counts(sub_con)

    add_constraint_scores(sub_con)

    add_select_score(sub_con)

    add_aggregate_overlapping_regions(sub_con)

    add_to_bed(sub_con)

    add_score_txt(sub_con)

    add_score_bed(sub_con)

    add_score_vcf(sub_con)

    ## Initiate argparse
    args = parser.parse_args(args)


    if args.check_config_path:

        print("Base Config Dir:\n", ", ".join(os.listdir(args.base_config)), "\n")
        print("ML Model Config Dir:\n", ", ".join(os.listdir(os.path.join(args.base_config,"ConSpliceML_Model"))), "\n")

        ## Check that the path exists
        assert os.path.exists(args.base_config), "!!ERROR!! The base config path does not exists. Bad path = '{}'\n".format(args.base_config)
        assert os.path.exists(args.config_path), "!!ERROR!! The config yaml path does not exists. Bad path = '{}'\n".format(args.config_path)
        assert os.path.exists(os.path.join(args.base_config,"ConSpliceML_Model")), "!!ERROR!! The ML Model config path does not exists. Bad path = '{}'\n".format(os.path.join(args.base_config,"ConSpliceML_Model"))
        assert os.path.exists(os.path.join(args.base_config,"ConSpliceML_Model/trained_ConSpliceML.rf")), "!!ERROR!! The ML Model .rf config path does not exists. Bad path = '{}'\n".format(os.path.join(args.base_config,"ConSpliceML_Model/trained_ConSpliceML.rf"))
        assert os.path.exists(os.path.join(args.base_config,"ConSpliceML_Model/training.yaml")), "!!ERROR!! The ML Model yaml config path does not exists. Bad path = '{}'\n".format(os.path.join(args.base_config,"ConSpliceML_Model/training.yaml"))

        ## Check that the file or dir is a file or dir
        assert os.path.isdir(args.base_config), "!!ERROR!! The base config is not a directory. Bad path = '{}'\n".format(args.base_config)
        assert os.path.isfile(args.config_path), "!!ERROR!! The config yaml is not a file. Bad path = '{}'\n".format(args.config_path)
        assert os.path.isdir(os.path.join(args.base_config,"ConSpliceML_Model")), "!!ERROR!! The ML Model config is not a directory. Bad path = '{}'\n".format(os.path.join(args.base_config,"ConSpliceML_Model"))
        assert os.path.isfile(os.path.join(args.base_config,"ConSpliceML_Model/trained_ConSpliceML.rf")), "!!ERROR!! The ML Model .rf config is not a file. Bad path = '{}'\n".format(os.path.join(args.base_config,"ConSpliceML_Model/trained_ConSpliceML.rf"))
        assert os.path.isfile(os.path.join(args.base_config,"ConSpliceML_Model/training.yaml")), "!!ERROR!! The ML Model yaml config is not a file. Bad path = '{}'\n".format(os.path.join(args.base_config,"ConSpliceML_Model/training.yaml"))

        ## Check that for non-empty files
        assert os.path.getsize(args.config_path) > 0, "!!ERROR!! The config yaml is empty. File Size = '{}'\n".format(os.path.getsize(args.config_path))
        assert os.path.getsize(os.path.join(args.base_config,"ConSpliceML_Model/trained_ConSpliceML.rf")), "!!ERROR!! The ML Model .rf config file is empty. File Size = '{}'\n".format(os.path.getsize(os.path.join(args.base_config,"ConSpliceML_Model/trained_ConSpliceML.rf")))
        assert os.path.isfile(os.path.join(args.base_config,"ConSpliceML_Model/training.yaml")), "!!ERROR!! The ML Model yaml config file is empty. File Size = '{}'\n".format(os.path.getsize(os.path.join(args.base_config,"ConSpliceML_Model/training.yaml")))

        print("\nAll config checks passed\n")

    elif args.command is None:
        
        print("\n!!ERROR!! Command not found. Please use the 'constraint', 'ML', '-h', or '--version' subcommand\n")

    else:

        try:
            args.func(parser, args)

        except AttributeError as e:

            if "'Namespace' object has no attribute 'func'" in str(e):
                print("\nUse '-h' to see the agrument options for this subcommand\n")

            else:
                print(e)

        

if __name__ == "__main__":
    sys.exit(main() or 0)
