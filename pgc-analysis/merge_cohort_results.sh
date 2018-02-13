#!/bin/bash

source params.sh

CHROM=$1

bcftools merge $(ls ${OUTDIR}/imputed/*bgs_imputeSTR_chr${CHROM}.vcf.gz) | \
    bgzip -c > ${OUTDIR}/imputed/PGC_imputeSTR_chr${CHROM}.vcf.gz
tabix -p vcf ${OUTDIR}/imputed/PGC_imputeSTR_chr${CHROM}.vcf.gz
