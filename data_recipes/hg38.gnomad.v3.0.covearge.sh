#!/bin/sh
set -eo pipefail -o nounset

genome=https://raw.githubusercontent.com/gogetdata/ggd-recipes/master/genomes/Homo_sapiens/hg38/hg38.genome

wget --quiet -O - https://storage.googleapis.com/gnomad-public/release/3.0/coverage/genomes/gnomad.genomes.r3.0.coverage.summary.tsv.bgz \
    | bgzip -d \
    | awk -v OFS="\t" '{ if (NR == 1 ) {print "#chrom","start","end",$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13} else if(NR > 1) {split($1,locus,":"); print locus[1],locus[2]-1,locus[2],$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}}' \
    | gsort /dev/stdin $genome \
    | bgzip -c > hg38.gnomad.v3.0.genome.coverage.summary.bed.gz

tabix hg38.gnomad.v3.0.genome.coverage.summary.bed.gz
