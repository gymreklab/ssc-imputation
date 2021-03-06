#!/bin/bash

source params.sh

#for chrom in $(seq 1 22)
#do
#chrom=20
#SNPSTR=${OUTDIR}/snpstrs/snpstrs_${chrom}.tab

SNPSTR=${OUTDIR}/na12878_snpstrs.bed

echo "chrom str snp snp10x snp1kg str10x str1kg"
while IFS='' read -r line || [[ -n "$line" ]]; do
    chrom=$(echo ${line} | cut -f 1 -d':')
    str=$(echo ${line} | cut -f 1 -d' ' | cut -f 2 -d':')
    snp=$(echo ${line} | cut -f 2 -d' ' | cut -f 2 -d':')
    if [[ x"$str" == x"locus1" ]]; then continue; fi;
    # Extract 10X SNP
    gt10x=$(bcftools query -r "chr${chrom}:${snp}-${snp}" -s NA12878_WGS_v2 -f "[%GT]\n" ${VCF10X} | head -n 1)
    if [[ x"$gt10x" == x"" ]]; then continue; fi;
    # Extract panel SNP
    gtpanel=$(bcftools query -r "${chrom}:${snp}-${snp}" -s NA12878 -f "[%GT]\n" ${KGDIR}/1kg.snp.str.chr${chrom}.vcf.gz | head -n 1)
    # Extract panel STR
    strpanel=$(bcftools query -r "${chrom}:${str}-${str}" -s NA12878 -f "%REF\t[%TGT]\n" ${KGDIR}/1kg.snp.str.chr${chrom}.vcf.gz | head -n 1)
    if [[ x"$strpanel" == x"" ]]; then continue; fi;
    reflen=$(echo $strpanel | cut -f 1 | awk '{print length($1)}')
    a1len=$(echo $strpanel | cut -f 2 | cut -f 1 -d'|' | awk '{print length($1)}')
    a2len=$(echo $strpanel | cut -f 2 | cut -f 2 -d'|' | awk '{print length($1)}')
    strpanel=$(($a1len-$reflen))"|"$(($a2len-$reflen))
    # Extract STR 10x
    str101=$(bcftools query -r "chr${chrom}:${str}-${str}" -f "[%GB]\n" ${OUTDIR}/NA12878_10x_phase1.vcf.gz | head -n 1)
    str102=$(bcftools query -r "chr${chrom}:${str}-${str}" -f "[%GB]\n" ${OUTDIR}/NA12878_10x_phase2.vcf.gz | head -n 1)
    str10x=${str101}"|"${str102}
    if [[ x"$str101" == x"" ]] || [[ x"$str102" == x"" ]]; then continue; fi;
    if [[ x"$str101" == x"." ]] || [[ x"$str102" == x"." ]]; then continue; fi;
    echo $chrom $str $snp $gt10x $gtpanel $str10x $strpanel
done < $SNPSTR


#done
