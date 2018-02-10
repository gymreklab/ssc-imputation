#!/bin/bash

#PBS -lwalltime=2:00:00 -lnodes=1:ppn=1
#PBS -d /home/gymrek/workspace/ssc-imputation/pgc-analysis

source params.sh

#CHROM=$1
CHROM=17

./pgc_regression.py \
    --strvcf ${OUTDIR}/imputed/PGC_imputeSTR_chr${CHROM}.vcf.gz \
    --sampledata ${OUTDIR}/pgc_regression_data.tab > pgc_regr_chr${CHROM}.tab