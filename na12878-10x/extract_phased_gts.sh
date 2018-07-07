#!/bin/bash

source params.sh

#for chrom in $(seq 1 22)
#do
chrom=20
SNPSTR=${OUTDIR}/snpstrs/snpstrs_${chrom}.tab

echo "chrom start end snp10x snp1kg str10x str1kg"
while IFS='' read -r line || [[ -n "$line" ]]; do
    str=$(echo ${line} | cut -f 1 -d' ' | cut -f 2 -d':')
    snp=$(echo ${line} | cut -f 2 -d' ' | cut -f 2 -d':')
    if [[ x"$str" == x"locus1" ]]; then continue; fi;
    # Extract 10X SNP
    gt10x=$(bcftools query -r "chr${chrom}:${snp}-${snp}" -s NA12878_WGS_v2 -f "[%GT]\n" ${VCF10X})
    if [[ x"$gt10x" == x"" ]]; then continue; fi;
    # Extract panel SNP
    gtpanel=$(bcftools query -r "${chrom}:${snp}-${snp}" -s NA12878 -f "[%GT]\n" ${KGDIR}/1kg.snp.str.chr${chrom}.vcf.gz)
    # Extract panel STR
    strpanel=$(bcftools query -r "${chrom}:${str}-${str}" -s NA12878 -f "%REF\t[%TGT]\n" ${KGDIR}/1kg.snp.str.chr${chrom}.vcf.gz)
    reflen=$(echo $strpanel | cut -f 1 | awk '{print length($1)}')
    a1len=$(echo $strpanel | cut -f 2 | cut -f 1 -d'|' | awk '{print length($1)}')
    a2len=$(echo $strpanel | cut -f 2 | cut -f 2 -d'|' | awk '{print length($1)}')
    strpanel=$(($a1len-$reflen))"|"$(($a2len-$reflen))
    # Extract STR 10x
    str10x=TODO
    echo $chrom $str $snp $gt10x $gtpanel $str10x $strpanel
done < $SNPSTR


#done
