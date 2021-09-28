#!/bin/sh
set -eo pipefail -o nounset

genome=https://raw.githubusercontent.com/gogetdata/ggd-recipes/master/genomes/Homo_sapiens/GRCh38/GRCh38.genome

## Get the chromomsome mapping file
chrom_mapping=$(ggd get-files grch38-chrom-mapping-ucsc2ensembl-ncbi-v1 --pattern "*.txt")

wget --quiet -O - https://storage.googleapis.com/gnomad-public/release/3.0/coverage/genomes/gnomad.genomes.r3.0.coverage.summary.tsv.bgz \
    | bgzip -d \
    | awk -v OFS="\t" '{ if (NR == 1 ) {print "#chrom","start","end",$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13} else if(NR > 1) {split($1,locus,":"); print locus[1],locus[2]-1,locus[2],$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}}' \
    | gsort --chromosomemappings $chrom_mapping /dev/stdin $genome \
    | bgzip -c > grch38.gnomad.v3.0.genome.coverage.summary.bed.gz

tabix grch38.gnomad.v3.0.genome.coverage.summary.bed.gz
