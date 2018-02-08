#!/bin/bash

source params.sh

COHORT=$1
CHROM=$2

java -jar ${BEAGLE} \
    gt=${OUTDIR}/imputed/chr${CHROM}/tmp/$(basename $COHORT).vcf.gz \
    ref=${OUTDIR}/refpanel/chr${CHROM}.str.snp.feb18.vcf.gz \
    out=${OUTDIR}/imputed/chr${CHROM}/$(basename $COHORT)_imputeSTR \
    niterations=0 # TODO use this?