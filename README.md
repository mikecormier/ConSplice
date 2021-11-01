# ConSplice

The Constrained Splicing (ConSplice) python module command line interface (CLI)


## Summary 

ConSplice consists of 

  1) a statistical model of splicing constarint accross protein-coding genes in the human genome

  2) a machine learning model that utlizes evidence of splicing from splicing predictions and the constraint of those events from ConSplice to improve the interpritation of alternative splicing variants.

ConSplice thus provides an improved method to identify deleterious alternative splicing variants that affect rare disease. 

Additionally, ConSplice allows one to identify splicing variants of interest outside of the exon-intron junction; those variants commonly ignored during clinical diangosis but vital for proper splicing.  


### Constarint 

**ConSplice** is a statistical model created to predict that genetic constraint against cryptic splicing in the Human Genome. 

ConSplice uses population levels of purifying selection across coding and non-coding regions of protein-coding genes from ostensibly healty inviduals in [gnomAD](https://gnomad.broadinstitute.org/), [Karczewski et al., Nature 2020](https://www.nature.com/articles/s41586-020-2308-7), to infer genetic constraint. It also uses per-nucleotide splicing predictions from [SpliceAI](https://github.com/Illumina/SpliceAI), [Jaganathan et al., Cell 2019](https://www.sciencedirect.com/science/article/pii/S0092867418316295?via%3Dihub), to identifying regions of genes likely important for alternative splicing. Combining signals of purifying selection from gnomAD and splicing from SpliceAI allowes use to infer the constraint against cryptic splicing. 

The ConSplice module:
  1) Creates a substitution rate matrix using gnomAD and SpliceAI to build an expectation of splicing mutation for regions of protein coding genes. 
  2) Calcualtes observed variation from gnomAD and expected variation from the substitution matrix for all `regions` of protein-coding genes. 
  3) Creates on Observed (O) to Expected (E) ratio, O/E, which can be used to infer the deviation from expected splicing mutation for a given ratio. A small O/E score indicates there are a fewer observed splicing variants compared to expected splicing variants in that region and it is constrained against cryptic splicing. A large O/E score indicates there are more observed compared to expected splicing variants, and it is thus more unconstrained. 
  4) Transformed O/E scores into percentiles ranges from 0.0 to 1.0. A percentile score of 0.0 suggests the region is completely unconstrained (tolerant) against aberrant splicing while a percentile score of 1.0 suggest the region is completely constrained (intolerant) against aberrant splicing.  


### Ensemble Machine Learning Approach

**ConSpliceML** is a machine learning approach using the ConSplice constraint score with per-base predictions of splicing from [SpliceAI](https://github.com/Illumina/SpliceAI) and [SQUIRLS](https://squirls.readthedocs.io/en/latest/) to identify potential pathogenic alternative splicing variants.

The intuition behind this model is based on improving the interpritation of pathogenic splicing using the splicing constraint metric. By combining a per-nucleotide prediction of alternative splicing with the 
pathogenic interpritation model of constraint, we can improve both the identification and interpritation of alternative splicing variants in terms of rare disease. 

Both SpliceAI and SQUIRLS are used to support alternative splicing predictions, as well as identify alternative splicing events missed by one of the tools. 


## ConSplice Module 





## Data used by ConSplice

ConSplice requires multiple sources of data including gnomAD, SpliceAI, GENCODE, and others. Data curation scripts called *Data Recipes* have been created for transperancy and reproducability. 

Data recipes can be found in the `data_recipes` directory of this repository. 


## Running ConSplice
