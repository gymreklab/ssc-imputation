#!/bin/bash

source params.sh

# Extract heterozygous phased SNPs in NA12878
bcftools query -f "%CHROM\t%POS\t[%GT]\n" $VCF10X | \
    grep -v "/" | grep -v "0|0" | grep -v "1|1" | \
    awk '{print $1 "\t" $2 "\t" $2+1}' | > ${OUTDIR}/na12878_het_snps.bed

# Intersect with SNPs in our panel
for chrom in $(seq 1 22)
do
    cat ${OUTDIR}/na12878_het_snps.bed | sed 's/chr//' | \
	intersectBed -b ${KGDIR}/1kg.snp.str.chr${chrom}.vcf.gz -a stdin -f 1 
done

#for chrom in $(seq 1 22)
#do
#done

# To get best tag SNP
#chrom=20
#    cat ${BESTSNPDIR}/chr${chrom}_snp_str_ld_50k.txt | \
#	datamash -H groupby 1 max 7 -f | cut -f 1,2,7 > \
#	${OUTDIR}/snpstrs/snpstrs_${chrom}.tab
