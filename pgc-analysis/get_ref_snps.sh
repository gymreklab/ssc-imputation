#!/bin/bash

# Get list of SNP positions and alleles for the reference panel

source params.sh

CHROM=$1

REFPANEL=${OUTDIR}/refpanel/chr${CHROM}.str.snp.feb18.vcf.gz

zcat ${REFPANEL} | grep -v "^#" | cut -f 3 | grep -v ":" > ${OUTDIR}/refpanel/pos/ssc_refpanel_snps_chr${CHROM}.txt
zcat ${REFPANEL} | grep -v "^#" | awk '($3!~/:/)' | cut -f 1-5 > ${OUTDIR}/refpanel/pos/ssc_refpanel_snps_alleles_chr${CHROM}.txt
