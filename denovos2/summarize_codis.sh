#!/bin/bash

source params.sh

echo "**** By length *****"
cat ${OUTDIR}/*_bylength.locus_summary.tab | grep -v chrom | intersectBed -a stdin -b CODIS.bed -wa -wb | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t"$9 "\t" $10 "\t" $11 "\t" $23 "\t" $24 "\t" $19}'

echo "**** By seq *****"
cat ${OUTDIR}/*_byseq.locus_summary.tab | grep -v chrom | intersectBed -a stdin -b CODIS.bed -wa -wb | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t"$9 "\t" $10 "\t" $11 "\t" $23 "\t" $24 "\t" $19}'
