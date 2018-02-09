#!/bin/bash

source params.sh

CHROM=$1

# Get files list
CFILES=""
for COHORT in $(cat ${COHORTFILES})
do
    CFILES="${CFILES} ${OUTDIR}/imputed/chr${CHROM}/$(basename $COHORT)_imputeSTR_chr${CHROM}.vcf.gz"
done

# Merge
bcftools merge $CFILES | bgzip -c > ${OUTDIR}/imputed/PGC_imputeSTR_chr${CHROM}.vcf.gz
tabix -p vcf ${OUTDIR}/imputed/PGC_imputeSTR_chr${CHROM}.vcf.gz
