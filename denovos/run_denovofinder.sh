#!/bin/bash

source params.sh

OUTDIR=${BASEOUTDIR}/denovofinder

for chrom in $(seq ${startchrom} ${endchrom})
do
    VCFFILE=/storage/s1saini/hipstr_genomewide/chr${chrom}/hipstr_calls_${chrom}.vcf.gz
    SNPVCF=/storage/mgymrek/ssc-denovos/phased-snps/shapeit.chr${chrom}.vcf.gz
    ${DENOVOFINDER} \
	--str-vcf ${VCFFILE} \
	--fam ${FAMFILE} \
	--snp-vcf ${SNPVCF} \
	--denovo-vcf ${OUTDIR}/denovofinder_chr${chrom}.vcf.gz \
	--window-size ${WINDOWSIZE} \
	--max-alleles ${MAXALLELES}
    tabix -p vcf ${OUTDIR}/denovofinder_chr${chrom}.vcf.gz
done
