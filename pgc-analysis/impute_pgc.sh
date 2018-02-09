#!/bin/bash

source params.sh

COHORT=$1
CHROM=$2

echo "Imputation with:"
echo "ref: ${OUTDIR}/refpanel/chr${CHROM}.str.snp.feb18.vcf.gz"
echo "target: ${OUTDIR}/imputed/chr${CHROM}/tmp/$(basename $COHORT).vcf.gz"

java -jar ${BEAGLE} \
    gt=${OUTDIR}/imputed/chr${CHROM}/tmp/$(basename $COHORT).vcf.gz \
    ref=${OUTDIR}/refpanel/pos/ssc_refpanel_snps_alleles_chr${CHROM}.vcf.gz \
    out=${OUTDIR}/imputed/chr${CHROM}/tmp/$(basename $COHORT)_imputeSTR_chr${CHROM} \
    niterations=0 # TODO use this?

tabix -p vcf ${OUTDIR}/imputed/chr${CHROM}/tmp/$(basename $COHORT)_imputeSTR_chr${CHROM}.vcf.gz

# Extract STR only to final file
intersectBed -a ${OUTDIR}/imputed/chr${CHROM}/tmp/$(basename $COHORT)_imputeSTR_chr${CHROM}.vcf.gz \
    -b ${HIPREF} -header | \
    grep "^#\|STR_" | \
    bgzip -c > ${OUTDIR}/imputed/chr${CHROM}/$(basename $COHORT)_imputeSTR_chr${CHROM}.vcf.gz
tabix -p vcf ${OUTDIR}/imputed/chr${CHROM}/$(basename $COHORT)_imputeSTR_chr${CHROM}.vcf.gz