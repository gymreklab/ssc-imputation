#!/bin/bash

source params.sh

zcat ${OUTDIR}/ssc_autosomal_perlocus_observed.bed.gz | \
    head -n 1 > ${OUTDIR}/ssc_perlocus_train_intergenic.bed
zcat ${OUTDIR}/ssc_autosomal_perlocus_observed.bed.gz | grep -v start | \
    intersectBed -a stdin -b ${TRANSCRIBED} -v >> \
    ${OUTDIR}/ssc_perlocus_train_intergenic.bed
gzip -f ${OUTDIR}/ssc_perlocus_train_intergenic.bed
