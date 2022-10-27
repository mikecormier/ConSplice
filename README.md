# ConSplice

[![python-version](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/consplice/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![latest-Release Badge](https://img.shields.io/github/v/release/mikecormier/ConSplice?label=latest%20release%2Fversion)


The Constrained Splicing (ConSplice) python module command line interface (CLI)


## Summary 

ConSplice consists of 

  1) a statistical model of splicing constraint across protein-coding genes in the human genome

  2) a machine learning model that utilizes evidence of splicing from splicing predictions and the constraint of those events from ConSplice to improve the interpretation of alternative splicing variants.

ConSplice thus provides an improved method to identify deleterious alternative splicing variants that affect rare disease. 

Additionally, ConSplice allows one to identify splicing variants of interest outside of the exon-intron junction; those variants commonly ignored during clinical diagnosis but vital for proper splicing.  


#### Constraint 

**ConSplice** is a statistical model created to predict that genetic constraint against cryptic splicing in the Human Genome. 

ConSplice uses population levels of purifying selection across coding and non-coding regions of protein-coding genes from ostensibly healthy individuals in [gnomAD](https://gnomad.broadinstitute.org/), [Karczewski et al., Nature 2020](https://www.nature.com/articles/s41586-020-2308-7), to infer genetic constraint. It also uses per-nucleotide splicing predictions from [SpliceAI](https://github.com/Illumina/SpliceAI), [Jaganathan et al., Cell 2019](https://www.sciencedirect.com/science/article/pii/S0092867418316295?via%3Dihub), to identifying regions of genes likely important for alternative splicing. Combining signals of purifying selection from gnomAD and splicing from SpliceAI allows use to infer the constraint against cryptic splicing across a gene. 

The ConSplice module:
  1) Creates a substitution probability matrix using gnomAD and SpliceAI to build an expectation of splicing mutation for regions of protein coding genes. 
  2) Calculates observed variation from gnomAD and expected variation from the substitution matrix for all `regions` of protein-coding genes. 
  3) Creates on Observed (O) to Expected (E) ratio, O/E, which can be used to infer the deviation from expected splicing mutation for a given ratio. A small O/E score indicates there are a fewer observed splicing variants compared to expected splicing variants in that region and it is constrained against cryptic splicing. A large O/E score indicates there are more observed compared to expected splicing variants, and it is thus more unconstrained. 
  4) Transformed O/E scores into percentiles, ranging from 0.0 to 1.0. A percentile score of 0.0 suggests the region is completely unconstrained (tolerant) against aberrant splicing while a percentile score of 1.0 suggest the region is completely constrained (intolerant) against aberrant splicing.  


#### Ensemble Machine Learning Approach

**ConSpliceML** is a machine learning approach using the ConSplice constraint score with per-base predictions of splicing from [SpliceAI](https://github.com/Illumina/SpliceAI) and [SQUIRLS](https://squirls.readthedocs.io/en/latest/) to identify potential pathogenic alternative splicing variants at a per-nucleotide level.

The intuition behind this model is based on improving the interpretation of pathogenic splicing using the splicing constraint metric. By combining a per-nucleotide prediction of alternative splicing with the 
pathogenic interpretation model of constraint, we can improve both the identification and interpretation of alternative splicing variants in terms of rare disease. 

Both SpliceAI and SQUIRLS are used to support alternative splicing predictions, as well as identify alternative splicing events missed by one of the tools. 

> **_NOTE:_** For now, the ConSpliceML model requires ConSplice, SpliceAI, and SQUIRLS scores as parameters to train and score the model. In later versions the number of features in the vector will be dynamic based on user input. 



## Using the ConSplice CLI 

### Data used by ConSplice

ConSplice requires multiple sources of data including gnomAD, SpliceAI, GENCODE, and others. Data curation scripts called *Data Recipes* have been created for transparency and reproducability. 

Data recipes can be found in the [data_recipes](https://github.com/mikecormier/ConSplice/tree/main/data_recipes) directory of this repository. 

> **_NOTE:_** The data requirements outlined in the `data_recipes` directory, or similar data, is required prior to running the ConSplice CLI. 


### Installing ConSplice

#### Conda/Mamba is required in order to install ConSplice requirements

We suggest using [mamba](https://github.com/mamba-org/mamba) over conda as it provides improvements to the conda infrastructure while still utilizing the conda framework. 

Mamba can be installed into the base environment of a current conda installation: `conda install -c conda-forge mamba` 

However, we suggest that **mambforge** be installed in place of conda.

Mambaforge will set up a normal conda environment with mamba ready to go. Therefore you can use both mamba and conda interchangeably. That means any conda command,
such as `install`, `search`, `uninstall`, `create`, `info`, etc. you can use with mamba. 

**You can install mambaforge [here](https://github.com/conda-forge/miniforge#mambaforge).** 

If you would still prefer to use conda, we suggest using [miniconda](https://conda.io/en/latest/miniconda.html). 


#### Conda/Mamba versions

The installation of the ConSplice CLI was recently tested with the latest versions of conda and mamba. 
  - conda v4.14.0 
  - mamba v0.25.0

We suggest using these or newer versions of conda and mamba.


#### Conda/Mamba Channels

You will need to add a few conda channels to your configurations before you can install ConSplice. 

Use the following commands to add the required channels. 

```
conda config --add channels defaults
conda config --add channels ggd-genomics
conda config --add channels bioconda
conda config --add channels conda-forge
```

#### Installing on different platforms

The ConSplice CLI is a `noarch` python package, which means it is architecture free. This means, in theory, that it can be installed on any system. 

However, we have seen a few problems with users trying to install ConSplice onto newer OSX systems using the **Apple M1 Silicon chips**. This is because some of the 
packages ConSplice depends on are not supported by the new ARM64 chip architecture. 

Luckily, the M1 chips can natively run the older x86 infrastructure used by the Intel chips. Additionally, conda/mamba provides a way to install x86 packages onto an ARM64 M1 Silicon Mac. You just need to tell conda which architecture to use. 

Specifically, you will add a conda environment variable to the installation instructions. This is how you tell conda to use the osx x86 architecture:
`CONDA_SUBDIR=osx-64`

In the following examples, we show how to install ConSplice onto your system, including onto an M1 Silicon chip Mac with an ARM64 architecture. 


#### Installation Examples:

**Install using mamba**:

(The example below will create a new conda/mamba environment with python3 and the ConSplice CLI already installed in it. It will then activate the environment so the ConSplice CLI can be used.)

```
## Create a ConSplice conda environment with python3 and consplice packages installed
mamba create --name ConSplice python=3 consplice 

## Activate the ConSplice conda environment
mamba activate ConSplice
```

(The example below will create a new conda environment, activate the environment, and install the ConSplice CLI into that environment)
```
## Create a ConSplice conda environment with python3 installed
mamba create --name ConSplice python=3  

## Activate the ConSplice conda environment
mamba activate ConSplice

## Install ConSplice while in the new conda environment
mamba install -c bioconda consplice
```


> **_NOTE:_** The example above creates and installs the ConSplice CLI into a new conda environment. You must activate the new conda environment using `conda activate <environment name>` to use ConSplice.  


**Install using mamba on an M1 Mac with an ARM64 architecture**:

(The example below will create a new conda/mamba environment with python3 and the ConSplice CLI already installed in it using the x86 architecture on an M1 Mac with an ARM64 architecture. It will then activate the environment so the ConSplice CLI can be used and set the conda config to use the x86 architecture.)

```
## Create a ConSplice conda environment with python3 and consplice packages installed
CONDA_SUBDIR=osx-64 mamba create --name ConSplice_x86 python=3 consplice 

## Activate the ConSplice conda environment
mamba activate ConSplice_x86

mamba config --env --set subdir osx-64
```

**Install using conda**:

To install with conda, you will use the exact same commands as shown above while replacing `mamba` with `conda`.


**Install from GitHub**:

(Although you can install from GitHub, we recommend you install ConSplice using conda/mamba)

```
## Clone the ConSplice GitHub repo
git clone https://github.com/mikecormier/ConSplice

cd ConSplice

## Create a ConSplice conda environment with python3 installed
mamba create --name ConSplice python=3  

## Activate the ConSplice conda environment
mamba activate ConSplice

## install the required packages
mamba install --file requirements.txt

## Install the ConSplice module
python setup.py install
```

**Install from GitHub on an M1 Mac with an ARM64 architecture**:

(Although you can install from GitHub, we recommend you install ConSplice using conda/mamba)

```
## Clone the ConSplice GitHub repo
git clone https://github.com/mikecormier/ConSplice

cd ConSplice

## Create a ConSplice conda environment with python3 installed
CONDA_SUBDIR=osx-64 mamba create --name ConSplice_x86 python=3  

## Activate the ConSplice conda environment
mamba activate ConSplice_x86

## install the required packages
CONDA_SUBDIR=osx-64 mamba install --file requirements.txt

## Install the ConSplice module
python setup.py install
```






### Running ConSplice

ConSplice has two subcommands:
- `constraint`: The module for ConSplice
- `ML`: The module for ConSpliceML 


Shown here is the top level help message
```
consplice -h


    =====================================================
    || ConSplice-CLI: The CLI to the ConSplice project ||
    =====================================================

    (1) Score/Annotate variants using the ConSpliceML model to identify potential pathogenic clinically relevant alternative splicing variants.
    (2) Train a Random Forest model for deleterious splicing prediction and interpretation utilizing the constrained splicing model.
    (3) Generate genic or regional splicing constraint profiles.

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Installed ConSplice CLI Version
  --config-path CONFIG_PATH
                        (constraint subcommand only) File path to the new ConSplice config yaml file. By default, ConSplice CLI will use the
                        default ConSplice config.yml file

Main Commands:
  {ML,constraint}
    ML                  ConSpliceML - Score variants with ConSpliceML or train a new ConSpliceML model
    constraint          Generate genic or regional splicing constraint profiles

```


#### constraint subcommand

This subcommand provides the following options
- `sub-matrix`:          Create a substitution matrix based on population variant frequency  
- `oe-counts`:           Using the substitution matrix, calculate observed and expected counts from population variation
- `calculate-oe`:        Calculate the O/E and Percentile scores 
- `select-score`:        Select scores and filter score file
- `agg-overlapping-reg`: Aggregate O/E and Percentile score when overlapping regions occur. (Only use if the step size is smaller than the window size when using the **oe-counts** command) 
- `to-bed`:              Convert filtered score file to bed format
- `score-txt`:           Add ConSplice scores to a txt file
- `score-bed`:           Add ConSplice scores to a bed file
- `score-vcf`:           Add ConSplice scores to a vcf file

A normal workflow for generating constraint profiles would 

  1) Create a substitution matrix 
  2) Generate observed and expected counts
  3) Calculate O/E and Percentile scores
  4) Extract the O/E and Percentiles scores and filter for relevant data 
  5) Convert filtered file to bed format 


```
consplice sub-matrix <arguments>
            |
            V
consplice oe-counts <arguments>
            |
            V
consplice calculate-oe <arguments>
            |
            V
consplice select-score <arguments>
            |
            V
consplice to-bed <arguments>
```

Once a bed file has been created, the consplice scores for a variants in a txt file, bed file, or vcf file can be added using the `score-txt`, `score-bed`, or `score-vcf` commands, respectively.

Shown here is the constraint level help message
```
    #######################
    # Splicing Constraint #
    #######################

    Module to generate genic or regional
    splicing constraint using patterns of
    purifying selection and evidence of
    alternative splicing

positional arguments:
  {sub-matrix,oe-counts,calculate-oe,select-score,agg-overlapping-reg,to-bed,score-txt,score-bed,score-vcf}
    sub-matrix          Build a substitution matrix using gnomAD variants and SpliceAI scores
    oe-counts           Calculate Observed and Expected splicing variant counts for genic or intragenic regions of genes
    calculate-oe        Calculate the O/E and Percentile constraint scores
    select-score        Select an O/E and matching Percentile score to filter on and remove all other non-essential columns after ConSplice scoring
    agg-overlapping-reg
                        Aggregate scores from overlapping regions where the step size of a region is less than the window size of the region
    to-bed              Convert the 1-based scored ConSplice txt file to a 0-based bed file
    score-txt           Add ConSplice scores to a tab-delimited txt variant file
    score-bed           Add ConSplice scores to a bed variant file
    score-vcf           Add ConSplice scores to a vcf file

optional arguments:
  -h, --help            show this help message and exit
```

#### ML subcommand

This subcommand provides the following options 
- `train`:     Train the ConSpliceML model
- `score-vcf`: Score a vcf file using a trained ConSpliceML model

ConSpliceML comes with a pre-trained ConSpliceML model when installed. Therefore, unless the need for training a new model arises, the most common command used 
with the ML subcommand is `score-vcf`. 

> **_NOTE:_** In order to score a vcf file with ConSpliceML, the vcf file must first have SpliceAI, SQURILS, and ConSplice annotations for each variant. ConSpliceML cannot score a vcf file without these annotations   


Shown here is the ML level help message
```
    ###############
    # ConSpliceML #
    ###############

    Module to:
     - Score variants to determine pathogenic alternative splicing variants using the ConSpliceML model.
     - Train a new ConSpliceML model.

positional arguments:
  {train,score-vcf}
    train            Train a Random Forest model using ConSplice.
    score-vcf        Score variants using ConSpliceML

optional arguments:
  -h, --help         show this help message and exit


```

> **_NOTE:_** `score-vcf` can be slow for a large vcf file. For your convince, a vcf file with pre-computed ConSpliceML scores for simulated SNVs in protein-coding genes can be found [here](). This mitigates the need to run `score-vcf`   


## Pre-computed ConSpliceML scores for SNVs in protein-coding genes

We provided a vcf file with pre-computed ConSpliceML scores for simulated SNVs in protein-coding genes. 

This file contains a SpliceAI score, a SQUIRLS score, a ConSplice Score, and a ConSpliceML score for each variant. SpliceAI, SQURILS, and ConSplice are included as a reference for the scores used to calculate the ConSpliceML score.

This vcf file can be found [here](https://home.chpc.utah.edu/~u1138933/ConSplice/scored_vcf/).


