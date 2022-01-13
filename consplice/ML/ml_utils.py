from __future__ import print_function

import io
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier

# ---------------------------------------------------------------------------------------------------------------------------------
## Global Variables
# ---------------------------------------------------------------------------------------------------------------------------------

patho_label = "patho_label"

trained_model = "trained_ConSpliceML.rf"

ml_info_file = "training.yaml"


# ---------------------------------------------------------------------------------------------------------------------------------
## Functions/Methods
# ---------------------------------------------------------------------------------------------------------------------------------


def ConSpliceML_train(
    training_df,
    feature_col_names=[],
    label_col="label",
    n_estimators=100,
    random_state=123,
):
    """
    ConSpliceML_train
    =================
    Method to train the ConSpliceML model using a Random Forest from sklearn.

    There is no limit to the number of features that can be used as long as those
    features are present in the training df.

    NOTE: This function will fail if there are missing/na values in the dataset
     for any of the features. The missing/na values need to be removed from the
     dataset prior to training the model.

    Parameters:
    -----------
    1) training_df: (Pandas DataFrame) A pandas dataframe where each row represents a training sample, and columns represent the features to train on.
                                        NOTE: features not trained on can be included in the df, but will be filtered out if they are not in the feature
                                        list.
                                        NOTE: The feature names must be in the columns of the df along with the label column.
    2) feature_col_names:       (list) A list of feature names to use to train the model. Every single name in this list should be in the df. The order
                                        of the names in the list will be used for the feature input to the model.
    3) label_col:                (str) The name of the label column in the df. Labels should be 1 for pathogenic and 0 for benign
    4) n_estimators:             (int) The number of estimators, decision tress, to build the Random Forest with. Default = 100
    5) random_state:             (int) The random state to build the model with. Used for reproducibility.

    Returns:
    ++++++++
    1) A trained sklearn Random Forest model.
    """

    ## RF models
    rf = RandomForestClassifier(random_state=random_state, n_estimators=n_estimators)

    rf.fit(training_df[feature_col_names], training_df[label_col].to_numpy())

    return rf


def ConSpliceML_score(rf, var_df, feature_col_names, new_col_name="ConSpliceML"):

    """
    ConSpliceML_score
    =================
    Score a set of variants using the ConSpliceML Random Forest trained model.

    Parameters:
    -----------
    1) rf:          (sklearn RF) The ConSpliceML trained RF model
    2) var_df:       (Pandas DF) A pandas Dataframe with variant info. Must contain the feature columns used to train the model
    3) feature_col_names: (list) A sorted list of features that match the scored ConSpliceML model to used to score the variant df.
    4) new_col_name:       (str) The name of the new column to add to the dataframe

    Returns:
    ++++++++
    1) (Pandas DataFrame) An updated pandas dataframe with the ConSpliceML score for the variants in the file
    2)              (str) The name of the ConSpliceML column added to the df.
    """

    var_df[new_col_name] = rf.predict_proba(var_df[feature_col_names])[:, 1]
    scored_df = var_df

    return (scored_df, new_col_name)


def get_alternative_gene_symbols(gene_info_file):
    """
    get_alternative_gene_symbols
    ============================
    This method is used to parse a file with accepted and alternative gene symbols and create a mapping between the pair.
     This function is set up to work with the HGNC Database mapping file. http://www.genenames.org/. The HGNC has put together
     A large data file that contains mapping between genes across the different data providers and each of the unique ids that
     each provider makes. Here, the function will use the accepted symbols and the alternative symbols from this file
     to create a mapping between the symbols. File url: ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt

    Parameters:
    ----------
    1) gene_info_file: (str) The file path to the HGNC symbol mapping file. NOTE: The file needs a "#" at the beginning of the 1st line in order to work.

    Returns:
    +++++++
    1) (dictionary) A dict where keys represent a single gene symbol, and values represent a set of synonymous symbols. The keys will include the accepted and
                    alternative gene symbols.

    NOTE: Any gene with no alternative symbols will be included as a key with the value as an empty set
    """

    ## Open the alt gene file
    try:
        gene_info_fh = (
            gzip.open(gene_info_file, "rt", encoding="utf-8")
            if gene_info_file.endswith(".gz")
            else io.open(gene_info_file, "rt", encoding="utf-8")
        )
    except IOError as e:
        print(
            "!!ERROR!! Unable to open the Alt Gene Symbols file: {}".format(
                gene_info_file
            )
        )
        print(str(e))
        sys.exit(1)

    ## Main dictionary for alternative symbols
    alternative_gene_symbols = defaultdict(set)

    ## Get an index from the header
    header_line = gene_info_fh.readline()
    assert (
        "#" in header_line
    ), "The alternative gene symbol file does not have a header staring with '#'. A header is required. Please add a header and try again"
    header_index = header_line.strip().replace("#", "").split("\t")

    ## Iterate over the file
    for line in gene_info_fh:

        ## Get a dictionary for the current line
        line_dict = dict(zip(header_index, line.strip().split("\t")))

        ## Get the official gene symbol and the synonymous symbols
        gene_symbol = line_dict["symbol"]

        ## Get that alternative gene symbols
        synonymous_symbols = (
            set(
                x for x in line_dict["alias_symbol"].strip().replace('"', "").split("|")
            )
            if line_dict["alias_symbol"]
            else set()
        )
        ## Add any previous symbols
        synonymous_symbols.update(
            set(x for x in line_dict["prev_symbol"].strip().replace('"', "").split("|"))
            if line_dict["prev_symbol"]
            else set()
        )
        ## Add current symbol
        synonymous_symbols.add(gene_symbol)

        ## Add to dict
        alternative_gene_symbols[gene_symbol].update(synonymous_symbols)

        ## iterate over each synonymous symbol
        for synonymous_symbol in synonymous_symbols:

            ## add synonymous symbol
            alternative_gene_symbols[synonymous_symbol].update(
                set(
                    [x for x in synonymous_symbols if x != synonymous_symbol]
                    + [gene_symbol]
                )
            )

    ## Close fh
    gene_info_fh.close()

    return alternative_gene_symbols
