#!/bin/sh
set -eo pipefail -o nounset

## Download the alt gene name file, add "#" to the header line
wget --quiet -O - ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt \
    | awk '{ if (NR == 1) print "#"$0; else print $0}' \
    > alt_gene_names.txt
