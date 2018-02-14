#!/bin/bash

set -e

source params.sh

COHORT=$1
CHROM=$2

java -jar ${BEAGLE} \
    gt=${OUTDIR}/imputed/tmp/$(basename $COHORT)_chr${CHROM}.vcf.gz \
    ref=${OUTDIR}/refpanel/pos/ssc_refpanel_snps_alleles_pgcregions_chr${CHROM}.vcf.gz \
    out=${OUTDIR}/imputed/tmp/$(basename $COHORT)_imputeSTR_chr${CHROM} \
    niterations=0 # TODO use this?

tabix -p vcf ${OUTDIR}/imputed/tmp/$(basename $COHORT)_imputeSTR_chr${CHROM}.vcf.gz

# Extract STR only to final file
intersectBed -a ${OUTDIR}/imputed/tmp/$(basename $COHORT)_imputeSTR_chr${CHROM}.vcf.gz \
    -b ${HIPREF} -header | \
    awk '(length($4)!=1)' | \
    bgzip -c > ${OUTDIR}/imputed/$(basename $COHORT)_imputeSTR_chr${CHROM}.vcf.gz
tabix -p vcf ${OUTDIR}/imputed/$(basename $COHORT)_imputeSTR_chr${CHROM}.vcf.gz