#!/bin/bash
# Example: ./get_pgc_vcf.sh /home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/scz/wave2/v1/cobg_dir_genome_wide/scz_zhh1_eur-qc.bgs

source params.sh

COHORT=$1 # Plink filename prefix
CHROM=$2

plink \
    --bfile ${COHORT} \
    --recode vcf bgz \
    --out ${OUTDIR}/imputed/tmp/$(basename $COHORT)_chr${CHROM} \
    --extract ${OUTDIR}/refpanel/pos/ssc_refpanel_snps_pgcregions_chr${CHROM}.txt \
    --real-ref-alleles \
    --a2-allele ${OUTDIR}/refpanel/pos/ssc_refpanel_snps_alleles_pgcregions_chr${CHROM}.txt 4 3 '#'

tabix -p vcf ${OUTDIR}/imputed/tmp/$(basename $COHORT)_chr${CHROM}.vcf.gz