#!/bin/bash

source params.sh

# Extract heterozygous phased SNPs in NA12878
bcftools query -f "%CHROM\t%POS\t[%GT]\t%FILTER\n" $VCF10X | grep PASS | \
    grep -v "/" | grep -v "0|0" | grep -v "1|1" | \
    awk '{print $1 "\t" $2-1 "\t" $2}'  > ${OUTDIR}/na12878_het_snps.bed

# Intersect with SNPs in our panel
for chrom in $(seq 1 22)
do
    cat ${OUTDIR}/na12878_het_snps.bed | sed 's/chr//' | \
	intersectBed -a ${KGDIR}/1kg.snp.str.chr${chrom}.vcf.gz -b stdin -f 1 -wa -wb | \
	awk '($2==$(NF))' | 
	awk '{print $(NF-2) "\t" $(NF-1) "\t" $NF}'
done > ${OUTDIR}/na12878_het_snps_inpanel.bed
bedtools sort -i ${OUTDIR}/na12878_het_snps_inpanel.bed > ${OUTDIR}/na12878_het_snps_inpanel_sorted.bed

# Get closest SNP to each STR
closestBed -a ${REGIONS2} -b ${OUTDIR}/na12878_het_snps_inpanel_sorted.bed -d -io | \
    awk '{print $1":"$2 "\t" $7":"$8 "\t" $10}' | grep -v "X:" | grep -v "Y:" > \
    ${OUTDIR}/na12878_snpstrs.bed
 
