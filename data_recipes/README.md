# ConSplice Data Recipes

These recipes represent the curation scripts to obtain the required data for running the ConSplice module.

Data not requiring a *recipe* because it is available from an external source will be indicated as such.


## ConSplice Data requirements

ConSplice requires the following data to calculate constraint. 

| Data                                            | Required for ConSplice   | Required for ConSpliceML | Recipe Name                                                                                                                                | Source           |
| -----                                           | ----------------------   | ------------------------ | -----------                                                                                                                                | ------           |
| gnomAD vcf                                      |        YES               |          NO              | [hg38.gnomad.genomes.v3.1.1.sites.sh](https://github.com/mikecormier/ConSplice/blob/main/data_recipes/hg38.gnomad.genomes.v3.1.1.sites.sh) | [See recipe below](https://github.com/mikecormier/ConSplice/tree/main/data_recipes#gnomad-vcf) |
| gnomAD Coverage                                 |        YES               |          NO              | [hg38.gnomad.v3.1.1.covearge.sh](https://github.com/mikecormier/ConSplice/blob/main/data_recipes/hg38.gnomad.v3.1.1.covearge.sh)           | [See recipe below](https://github.com/mikecormier/ConSplice/tree/main/data_recipes#gnomad-coverage) |
| Alternative Gene Names                          |        YES               |          YES             | [alt_gene_names.sh](https://github.com/mikecormier/ConSplice/blob/main/data_recipes/alt_gene_names.sh) | [See recipe below](https://github.com/mikecormier/ConSplice/tree/main/data_recipes#alternativesynonymous-gene-names) |
| SpliceAI                                        |        YES               |  YES - score required    | NA | [See recipe below](https://github.com/mikecormier/ConSplice/tree/main/data_recipes#spliceai) |
| SQUIRLS                                         |        NO                |  YES - score required    | NA | [See recipe below](https://github.com/mikecormier/ConSplice/tree/main/data_recipes#squirls) |
| grch38-reference-genome-gencode-v1              |        YES               |          NO              | NA | [GGD](https://gogetdata.github.io/recipes/genomics/Homo_sapiens/GRCh38/grch38-reference-genome-gencode-v1/README.html) |
| grch38-canonical-transcript-features-gencode-v1 |        YES               |          NO              | NA | [GGD](https://gogetdata.github.io/recipes/genomics/Homo_sapiens/GRCh38/grch38-canonical-transcript-features-gencode-v1/README.html) |
| grch38-segmental-dups-ucsc-v1                   |        YES               |          NO              | NA | [GGD](https://gogetdata.github.io/recipes/genomics/Homo_sapiens/GRCh38/grch38-segmental-dups-ucsc-v1/README.html) |
| grch38-self-chain-ucsc-v1                       |        YES               |          NO              | NA | [GGD](https://gogetdata.github.io/recipes/genomics/Homo_sapiens/GRCh38/grch38-self-chain-ucsc-v1/README.html) |
| Pathogenic truth set                            | YES + filtering required |          YES             | NA | [see recipe below](https://github.com/mikecormier/ConSplice/tree/main/data_recipes#truth-sets) |
| Benign truth set                                | YES + filtering required |          YES             | NA | [see recipe below](https://github.com/mikecormier/ConSplice/tree/main/data_recipes#truth-sets) |



`YES - score required` indicates that the score is needed but that module does not use the score from the primary data source. For example, ConSpliceML requires SpliceAI scores to be annotated into the vcf file but does not use the scores from the primary SpliceAI file. 
`YES - filtering rqeuired` indicates that the data is used directly by the model and should be used to filter other data that is used by the model. For example, pathogenic and benign variants used to test the model should be used to filter variants out of gnomAD for training the model. 

## Data Recipe Description:


### gnomAD VCF:

> **_NOTE:_** The recipe for gnomAD v3.0 variants is not provided here, but the following script could be modified with urls for the v3.0 version instead of the v3.1.1 version provided below.  


#### hg38.gnomad.genomes.v3.1.1.sites.sh 

A recipe to download the per chromosome vcf files for gnomAD v3.1.1 with 76,156 whole genomes. Once downloaded, these files will be combined, converted to *bcf* format and index using a *csi* index. It will clean up any intermediate files once the final gnomAD bcf file is created. 

The *hg38* designation is used here since these files contain a *chr* prefix. 

**As a note, this process requires a large amount of storage.** 

The final data file will be named *hg38.gnomad.genomes.v3.1.1.sites.bcf* with a csi index *hg38.gnomad.genomes.v3.1.1.sites.bcf.csi*

the final data file will be ~1T in size

This script produces a **bcf** for processing. ConSplice can processes bcf files much faster than vcf files.


#### Filtered hg38.gnomad.genomes.v3.1.1.sites.bcf

We further filter the *hg38.gnomad.genomes.v3.1.1.sites.bcf* file using all variants in the truth sets. (See the explanation of truth sets below) 

We filter these variants out of gnomAD prior to running ConSplice to remove the set of variants we will test both ConSplice and ConSpliceML on, and thus,
remove the bias of training and testing on the same set of data. 

To filter this data out, we use the `isec` function of `bcftools` and select the resulting file the contains the private gnomAD variants. 

```
bcftools isec -p <output_dir> -O b <combined-truth-set-vcf> hg38.gnomad.genomes.v3.1.1.sites.bcf

```

The final gnomAD file with the truth set variants removed is the file used for ConSplice.



### gnomAD Coverage: 

A few of these recipes require other GGD recipes to be download prior to running the recipe script. To run this script install the required GGD data packages from the list above and [grch38-chrom-mapping-ucsc2ensembl-ncbi-v1](https://gogetdata.github.io/recipes/genomics/Homo_sapiens/GRCh38/grch38-chrom-mapping-ucsc2ensembl-ncbi-v1/README.html) 

#### grch38.gnomad.v3.0.covearge.sh and hg38.gnomad.v3.0.covearge.sh

> **_NOTE:_**  hg38.gnomad.v3.0.covearge.sh uses the *chr* prefix while grch38.gnomad.v3.0.covearge.sh does not use the *chr* prefix.

A recipe to download and curate the gnomAD **v3.0** GRCh38/hg38 coverage file used by ConSplice to identify gnomAD with adequate sequencing coverage.

This coverage file will be downloaded from the gnomAD cloud storage, processed and transformed into bed format, bgzip, and tabixed.  

Data file size for the output file *GRCh38.gnomad.v3.0.genome.coverage.summary.bed.gz* or *hg38.gnomad.v3.0.genome.coverage.summary.bed.gz* is ~70G after processing 


#### grch38.gnomad.v3.1.1.covearge.sh and hg38.gnomad.v3.1.1.covearge.sh

> **_NOTE:_**  hg38.gnomad.v3.1.1.covearge.sh uses the *chr* prefix while grch38.gnomad.v3.1.1.covearge.sh does not use the *chr* prefix.

A recipe to download and curate the gnomAD **v3.1.1** GRCh38/hg38 coverage file used by ConSplice to identify gnomAD with adequate sequencing coverage.

This coverage file will be downloaded from the gnomAD cloud storage, processed and transformed into bed format, bgzip, and tabixed.  

Data file size for the output file *GRCh38.gnomad.v3.1.1.genome.coverage.summary.bed.gz* or *hg38.gnomad.v3.1.1.genome.coverage.summary.bed.gz* is ~80G after processing 



> **_NOTE:_** for ConSplice, we use the hg38.gnomad.v3.#.genome.coverage.summary.bed.gz file with the *chr* prefix which matches the gnomAD v# vcf file 



### Alternative/Synonymous Gene Names

There may be some inconsistencies between the gene name in gnomAD, GENCODE, and SpliceAI. Therefore, we use a file that provides alternative/synonymous gene symbols for the same gene to check for matching genes.

#### alt_gene_names.sh

This recipe downloads and formats the alternative gene symbol file for use with ConSplice. The final data file will be named *alt_gene_names.txt*




## External data files

### GGD Data Recipes

All GGD data packages used by ConSplice are available through [GGD](https://gogetdata.github.io/). The GGD cli should be used to download and use these data files.


#### GGD recipe list:

  - `grch38-chrom-mapping-ucsc2ensembl-ncbi-v1`  (This recipe is used to support the processing of other recipes listed above)

  - `grch38-reference-genome-gencode-v1` 
    
    This reference genome comes bgzipped. To improve speed, we recommend decompressing the file prior to using it with ConSplice

  - `grch38-canonical-transcript-features-gencode-v1`

  - `grch38-segmental-dups-ucsc-v1`

  - `grch38-self-chain-ucsc-v1`



### SpliceAI 

ConSplice and ConSpliceML use the GRCh38 **RAW** SpliceAI SNV vcf file (not the masked file) named *spliceai_scores.raw.snv.hg38.vcf.gz*. 

You can download the SpliceAI variant file from illumina's [basespace](https://basespace.illumina.com/). 

Further details can be found on the [SpliceAI](https://github.com/Illumina/SpliceAI) github page. 

> **_NOTE:_** ConSpliceML does not use these scores directly from the *spliceai_scores.raw.snv.hg38.vcf.gz* but rather uses the scores that are annotated in a vcf file for each variant. To use ConSpliceML, SpliceAI scores needed to be added to the variant file prior to running ConSpliceML 


### SQUIRLS 

ConSpliceML uses the GRCh38 SQUIRLS SNV scores. SQUIRLS does not provided a pre-computed scores. Therefore, SQUIRLS needs to be run on the variant file prior to running ConSpliceML.

The SQUIRLS CLI can be found [here](https://squirls.readthedocs.io/en/latest/)

ConSpliceML can be run once SQUIRLS has been downloaded and SQUIRLS has scored a vcf file. 


### Truth Sets

Truth sets used to test ConSplice and ConSpliceML which should be filtered out of the gnomAD vcf

#### Pathogenic truth set

Deleterious splicing variants with a **DM** tag annotated in the HGMD database were used a pathogenic variants in the truth set. Due to licensing agreements with HGMD, we are not allowed to share these variants publicly. 

If you have an HGMD license you can obtain the data using the following steps:

```
  # 1) Download the HGMD SQL Dump. (The current model uses the 2021 q1 SQL dump)

  # 2) Use the following SQL command to extract deleterious splice altering variants 

    SELECT * FROM splice INNER JOIN hgmd_hg38_vcf ON acc_num = hgmd_hg38_vcf.id
```

#### Benign truth set

The benign truth set used in the current model can be found in Supplemental Table S5 of the ConSplice paper.

