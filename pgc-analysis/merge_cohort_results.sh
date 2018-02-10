#!/bin/bash

source params.sh

CHROM=$1

bcftools merge $(ls ${OUTDIR}/imputed/chr${CHROM}/*.vcf.gz | grep -v scz_clm2_eur-qc | grep -v scz_umeb_eur-qc ) | bgzip -c > ${OUTDIR}/imputed/PGC_imputeSTR_chr${CHROM}.vcf.gz
tabix -p vcf ${OUTDIR}/imputed/PGC_imputeSTR_chr${CHROM}.vcf.gz
