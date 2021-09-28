# Additional Data Recipes

These recipes represent the additional data recipe for data curation when such data recipes don't exist in GGD and the data is not explicitly available from another source


## Data Recipes:


### gnomAD VCF:

> **_NOTE:_** The recipe for gnomAD v3.0 variants is not provided here, but the following script could be modified with urls for the v3.0 version instead of the v3.1.1 version provided below.  


#### hg38.gnomad.genomes.v3.1.1.sites.sh 

A recipe to download the per chromsome vcf files for gnomAD v3.1.1 with 76,156 whole genomes. Once downloaded, these files will be combined, converted to *bcf* format and index using a *csi* index. It will clean up any intermediate files once the final gnomAD bcf file is created. 

The *hg38* designation is used here since these files contain a *chr* prefix. 

As a note, this process requires a lare amount of storage. 

The final data file will be named *hg38.gnomad.genomes.v3.1.1.sites.bcf* with a csi index *hg38.gnomad.genomes.v3.1.1.sites.bcf.csi*

the final data file will be ~1T in size


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

A few of these recipes require other GGD recipes to be download prior to runing the recipe script. See the GGD recipe list below.

#### grch38.gnomad.v3.0.covearge.sh and hg38.gnomad.v3.0.covearge.sh

> **_NOTE:_**  hg38.gnomad.v3.0.covearge.sh uses the *chr* prefix while grch38.gnomad.v3.0.covearge.sh does not use the *chr* prefix.

A recipe to download and currate the gnomAD **v3.0** GRCh38/hg38 coverage file used by ConSplice to identify gnomAD with adequate sequencing coverage.

This coverage file will be downloaded from the gnomAD cloud storage, processed and transformed into bed format, bgzip, and tabixed.  

Data file size for the output file *GRCh38.gnomad.v3.0.genome.coverage.summary.bed.gz* or *hg38.gnomad.v3.0.genome.coverage.summary.bed.gz* is ~70G after processing 


#### grch38.gnomad.v3.1.1.covearge.sh and hg38.gnomad.v3.1.1.covearge.sh

> **_NOTE:_**  hg38.gnomad.v3.1.1.covearge.sh uses the *chr* prefix while grch38.gnomad.v3.1.1.covearge.sh does not use the *chr* prefix.

A recipe to download and currate the gnomAD **v3.1.1** GRCh38/hg38 coverage file used by ConSplice to identify gnomAD with adequate sequencing coverage.

This coverage file will be downloaded from the gnomAD cloud storage, processed and transformed into bed format, bgzip, and tabixed.  

Data file size for the output file *GRCh38.gnomad.v3.1.1.genome.coverage.summary.bed.gz* or *hg38.gnomad.v3.1.1.genome.coverage.summary.bed.gz* is ~80G after processing 



> **_NOTE:_** for ConSplice, we use the hg38.gnomad.v3.#.genome.coverage.summary.bed.gz file with the *chr* prefix which matches the gnomAD v# vcf file 



### Alternative/Synonymous Gene Names

There may be some inconcistencies between the gene name in gnomAD, GENCODE, and SpliceAI. Therefore, we use a file that provides alternative/synonymous gene symbols for the same gene to check for matching genes.

#### alt_gene_names.sh

This recipe downloads and formats the altnerative gene symbol file for use with ConSplice. The final data file will be named *alt_gene_names.txt*



## GGD Data Recipes

All other recipes used by ConSplice are available throuhg [GGD](https://gogetdata.github.io/). The GGD cli should be used to donwload and use these data files.


### GGD recipe list:

  - `grch38-chrom-mapping-ucsc2ensembl-ncbi-v1`  (This recipe is used to support the processing of other recipes listed above)

  - `grch38-reference-genome-gencode-v1` 
    
    This reference genome comes bgzipped. To improve speed, we recommend decompressing the file prior to using it with ConSplice

  - `grch38-canonical-transcript-features-gencode-v1`

  - `grch38-segmental-dups-ucsc-v1`

  - `grch38-self-chain-ucsc-v1`



## SpliceAI 

ConSplice uses the GRCh38 **RAW** SpliceAI SNV vcf file (not the masked file) named *spliceai_scores.raw.snv.hg38.vcf.gz*. 

You can download the SpliceAI variant file from illumina's [basespace](https://basespace.illumina.com/). 

Further details can be found on the [SpliceAI](https://github.com/Illumina/SpliceAI) github page. 



## Truth Sets

[how to get the data, if there is a paper associated with the data, how to filter the data, etc.]


HGMD 

Validated GTEx from SpliceAI 

CEPH 

deCODE 

