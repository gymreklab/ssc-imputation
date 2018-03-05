#!/bin/bash
#SBATCH -A csd568
#SBATCH --get-user-env 
#SBATCH --job-name Create_Filters_Denovos
#SBATCH -p shared
#SBATCH -t 1000
#SBATCH -o /oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered/logs/Create_Filter_SSC_Phase1_Denovos.log
#SBATCH -e /oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered/logs/Create_Filter_SSC_Phase1_Denovos.log 


source params.sh

tmpdir=$(mktemp -d)

########################################
# Which loci to filter

# Clean loci filters
rm -f ${FILTERDIR}/locus_filter*.bed

# SSC VCF filters
for chrom in `seq 1 22`
do
    bcftools query -f "%CHROM\t%POS\t%INFO/END\t%FILTER\n" ${VCFDIR}/hipstr.chr${chrom}.allfilters.vcf.gz | \
    awk '($4!="PASS" && $4!="Het") {print $1 "\t" $2 "\t" $3}' > ${FILTERDIR}/locus_filter_ssc_filters_chr${chrom}.bed
done

# Combine all loci filters
cat ${FILTERDIR}/locus_filter*.bed | sort -k1,1 -k2,2n | uniq > \
    ${FILTERDIR}/denovo_locus_filters.bed

########################################
# Which families to filter

# Clean family filters
rm -f ${FILTERDIR}/family_filter*.txt

# This family had way too many denovos called in the affected child
# Use sample IDs instead of family IDs, which are just numbers
echo "SSC02229" > ${FILTERDIR}/family_filters_outlier.txt
echo "SSC02241" >> ${FILTERDIR}/family_filters_outlier.txt

# Combine all family filters
cat ${FILTERDIR}/*.txt | sort | uniq > ${FILTERDIR}/denovo_family_filters.txt


