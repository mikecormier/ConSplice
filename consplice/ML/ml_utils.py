from __future__ import print_function
import sys
import os
import io
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier


#---------------------------------------------------------------------------------------------------------------------------------
## Global Variables
#---------------------------------------------------------------------------------------------------------------------------------

patho_label = "patho_label"

trained_model = "trained_ConSpliceML.rf"

ml_info_file = "training.yaml"


#---------------------------------------------------------------------------------------------------------------------------------
## Functions/Methods
#---------------------------------------------------------------------------------------------------------------------------------


def ConSpliceML_train(training_df,
                      feature_col_names = [],
                      label_col = "label",
                      n_estimators = 100,
                      random_state = 123):
    """
    ConSpliceML_train
    =================
    Method to train the ConSpliceML model using a Random Forest from sklearn.

    There is no limit to the number of features that can be used as long as those 
    features are present in the training df. 

    NOTE: This function will fail if there are missing/na values in the dataset 
     for any of the features. The missing/na values need to be removed from the 
     dataset prior to trainig the model.

    Parameters:
    -----------
    1) training_df: (Pandas DataFrame) A pandas dataframe where each row represents a training sample, and columns represent the features to train on.
                                        NOTE: features not trained on can be included in the df, but will be filtered out if they are not in the feature 
                                        list. 
                                        NOTE: The feature names must be in the columns of the df along with the label column.
    2) feature_col_names:       (list) A list of feature names to use to train the model. Every single name in this list should be in the df. The order 
                                        of the names in the list will be used for the feature input to the model.
    3) label_col:                (str) The name of the label column in the df. Labels should be 1 for pathogenic and 0 for bening
    4) n_estimators:             (int) The number of estimators, decision tress, to build the Random Forest with. Default = 100
    5) random_state:             (int) The random state to build the model with. Used for reproducibility. 

    Returns:
    ++++++++
    1) A trained sklearn Random Forest model. 
    """

    ## RF models
    rf = RandomForestClassifier(random_state = random_state, n_estimators = n_estimators)

    rf.fit(training_df[feature_col_names],training_df[label_col].to_numpy())
    
    return(rf)
