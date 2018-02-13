#!/bin/bash

# Get list of SNP positions and alleles for the reference panel

source params.sh
CHROM=$1

#REFPANEL=${OUTDIR}/refpanel/chr${CHROM}.str.snp.feb18.vcf.gz

cat ${REGIONS} | awk -v"window=$WINDOW" '{print $1 "\t" $2-window "\t" $3+window}' | \
    intersectBed -a ${REFPANEL} -b stdin | grep -v "^#" | grep -w "^${CHROM}" | cut -f 3 | grep -v ":" \
    > ${OUTDIR}/refpanel/pos/ssc_refpanel_snps_pgcregions_chr${CHROM}.txt

cat ${REGIONS} | awk -v"window=$WINDOW" '{print $1 "\t" $2-window "\t" $3+window}' | \
    intersectBed -a ${REFPANEL} -b stdin | grep -v "^#" | grep -w "^${CHROM}" | \
    awk '($3!~/:/)' | cut -f 1-5 > ${OUTDIR}/refpanel/pos/ssc_refpanel_snps_alleles_pgcregions_chr${CHROM}.txt

tabix -h ${REFPANEL} ${CHROM} | bgzip -c > ${OUTDIR}/refpanel/pos/ssc_refpanel_snps_alleles_pgcregions_chr${CHROM}.vcf.gz
tabix -p vcf ${OUTDIR}/refpanel/pos/ssc_refpanel_snps_alleles_pgcregions_chr${CHROM}.vcf.gz

#cat ${REGIONS} | awk -v"window=$WINDOW" '{print $1 "\t" $2-window "\t" $3+window}' | \
#    intersectBed -a ${REFPANEL} -b stdin -header | bgzip -c > ${OUTDIR}/refpanel/pos/ssc_refpanel_snps_alleles_chr${CHROM}.vcf.gz
#tabix -p vcf ${OUTDIR}/refpanel/pos/ssc_refpanel_snps_alleles_chr${CHROM}.vcf.gz