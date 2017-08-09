#!/bin/bash

source params.sh

OUTDIR=${BASEOUTDIR}/denovofinder_trio

for chrom in $(seq ${startchrom} ${endchrom})
do
    VCFFILE=/storage/s1saini/hipstr_genomewide/chr${chrom}/hipstr_calls_${chrom}.vcf.gz
    ${DENOVOFINDER} \
	--str-vcf ${VCFFILE} \
	--fam ${FAMFILE} \
	--denovo-vcf ${OUTDIR}/denovofinder_chr${chrom}.vcf.gz \
	--max-alleles ${MAXALLELES}
    tabix -p vcf ${OUTDIR}/denovofinder_chr${chrom}.vcf.gz
done
