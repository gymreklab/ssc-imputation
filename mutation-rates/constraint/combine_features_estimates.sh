#!/bin/bash

source params.sh

tmpdir=$(mktemp -d)

zcat ${OUTDIR}/GRCh37.hipstr_reference_sorted_properties.tab | grep -v chrom | \
    intersectBed -a ${ESTDIR}/ssc_mutea_auto_scaled.bed.gz -b stdin -wa -wb | \
    cut -f 10-12 --complement > ${tmpdir}/tmp
echo "chrom,start,end,ml_mu,ml_beta,ml_p,ml_mu_stderr,numsamples,center,motif,length,uninterrupted_length" | \
    sed 's/,/\t/g' > ${OUTDIR}/ssc_autosomal_perlocus_observed.bed
cat ${tmpdir}/tmp >> ${OUTDIR}/ssc_autosomal_perlocus_observed.bed
gzip -f ${OUTDIR}/ssc_autosomal_perlocus_observed.bed
