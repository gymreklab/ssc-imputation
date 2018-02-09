#!/bin/bash

source params.sh

CHROM=$1

./pgc_regression.py \
    --strvcf ${OUTDIR}/imputed/PGC_imputeSTR_chr${CHROM}.vcf.gz \
    --sampledata ${OUTDIR}/pgc_regression_data.tab