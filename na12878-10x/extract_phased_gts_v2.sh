#!/bin/bash

source params.sh

SNPSTR=${OUTDIR}/na12878_snpstrs.bed

cat ${SNPSTR} | cut -f 2 | awk -F":" '{print "chr"$1 "\t" $2 "\t" $2+1}' | \
    grep -v "\." | uniq > ${OUTDIR}/tmp/snps_chr.bed
cat ${SNPSTR} | cut -f 2 | awk -F":" '{print $1 "\t" $2 "\t" $2+1}' | \
    grep -v "\." | uniq > ${OUTDIR}/tmp/snps.bed
cat ${SNPSTR} | cut -f 1 | awk -F":" '{print "chr"$1 "\t" $2 "\t" $2+1}' | \
    grep -v "\." | uniq > ${OUTDIR}/tmp/strs_chr.bed
cat ${SNPSTR} | cut -f 1 | awk -F":" '{print $1 "\t" $2 "\t" $2+1}' | \
    grep -v "\." | uniq > ${OUTDIR}/tmp/strs.bed

######### Extract STRs #############
# From 10X
intersectBed -a ${OUTDIR}/NA12878_10x_phase1.vcf.gz -b ${OUTDIR}/tmp/strs_chr.bed -header | \
    bcftools query -f "%CHROM\t%POS\t[%GB\n]" | awk '{print $1 "\t" $2 "\t" $2+1 "\t" $3}'  \
    > ${OUTDIR}/tmp/str10x_1.tab
intersectBed -a ${OUTDIR}/NA12878_10x_phase2.vcf.gz -b ${OUTDIR}/tmp/strs_chr.bed -header | \
    bcftools query -f "%CHROM\t%POS\t[%GB\n]" | awk '{print $1 "\t" $2 "\t" $2+1 "\t" $3}'  \
    > ${OUTDIR}/tmp/str10x_2.tab
intersectBed -a ${OUTDIR}/tmp/str10x_1.tab -b ${OUTDIR}/tmp/str10x_2.tab -wa -wb | \
    awk '($2==$6)' | awk '{print $1 "\t" $2 "\t" $3 "\t" $4"|"$8}' | \
    sed 's/chr//' | grep -v "\."  > ${OUTDIR}/tmp/str10x.tab

# From panel
for chrom in $(seq 1 22)
do
    intersectBed -a ${KGDIR}/1kg.snp.str.chr${chrom}.vcf.gz -b ${OUTDIR}/tmp/strs.bed -header | \
	bcftools view --max-alleles 20 | \
	bcftools query -s"NA12878" -f "%CHROM\t%POS\t%REF\t[%TGT\n]" | sed 's/|/\t/' | \
	awk '{print $1 "\t" $2 "\t" $2+1 "\t" (length($4)-length($3))"|"(length($5)-length($3))}'
done > ${OUTDIR}/tmp/strpanel.tab

# Merge
intersectBed -a $OUTDIR/tmp/strpanel.tab -b $OUTDIR/tmp/str10x.tab -wa -wb | \
    awk '($2==$6)' | cut -f 1-4,8 > ${OUTDIR}/tmp/merged_strs.bed

######### Extract SNPs #############
# Extract SNPs 10x
intersectBed -a ${VCF10X} -b ${OUTDIR}/tmp/snps_chr.bed -header -f 1 | \
    bcftools query -f"%CHROM\t%POS\t[%GT]\n" | sed 's/chr//' | \
    awk '{print $1 "\t" $2 "\t" $2+1 "\t" $3}' > ${OUTDIR}/tmp/gt10x.tab

# Extract SNPs panel
for chrom in $(seq 1 22)
do
    vcf=${KGDIR}/1kg.snp.str.chr${chrom}.vcf.gz
    intersectBed -a ${vcf} -b ${OUTDIR}/tmp/snps.bed -header -f 1| \
	bcftools query -s "NA12878" -f"%CHROM\t%POS\t[%GT]\n" | \
	awk '{print $1 "\t" $2 "\t" $2+1 "\t" $3}' 
done > ${OUTDIR}/tmp/gtpanel.tab

# Merge
intersectBed -a $OUTDIR/tmp/gtpanel.tab -b $OUTDIR/tmp/gt10x.tab -wa -wb | \
    awk '($2==$6)' | cut -f 1-4,8 > ${OUTDIR}/tmp/merged_snps.bed

# Merge everything to single file
echo "chrom,strpos,snppos,strpanel,str10x,snppanel,snp10x" | sed 's/,/\t/g' > ${OUTDIR}/na12878_panel_vs_10x.tab
cat $SNPSTR | sed 's/ /\t/g' | sed 's/:/\t/g' | \
    awk '{print $1 "\t" $2 "\t" $2+1"\t" $4}' | \
    intersectBed -a ${OUTDIR}/tmp/merged_strs.bed -b stdin -wa -wb | \
    awk '($2==$7)' | awk '{print $1 "\t" $9+1 "\t" $9+2 "\t" $2 "\t" $4 "\t" $5}' | \
    intersectBed -a stdin -b ${OUTDIR}/tmp/merged_snps.bed -wa -wb | \
    awk '{print $1 "\t" $4 "\t" $2-1 "\t" $5 "\t" $6 "\t" $10 "\t" $11}' \
    >> ${OUTDIR}/na12878_panel_vs_10x.tab 
