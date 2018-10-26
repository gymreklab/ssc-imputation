#!/bin/bash

COHORTFILES=/home/gymrek/pgc_imputation/pgc_eur_mergelist.txt
IMPUTEDVCF=/home/gymrek/ssaini/per_locus/11.57523883/PGC.11.57523883.imputed.noAF.vcf.gz
SCRATCHDIR=tmp/
CHROM=11
POS=57523883
COVARFILE=/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/scz/wave2/v1/prune.bfile.cobg.PGC_SCZ49.sh2.menv.mds_cov
COVARCOLS=C1,C2,C3,C4,C5,C6,C7,C9,C15,C18 # based on recommendation from S. Ripke
OUTDIR=.
SNPID=rs9420

# Get list of SNP columns
snpcols=$(zcat $IMPUTEDVCF | \
    grep -w $SNPID | \
    grep -v "^\#" | \
    awk '{print "SNP_"$2"_"$3","}' | \
    tr -d '\n' | sed 's/,$//')

# Extract all SNP genotypes and sample names to a file
zcat $IMPUTEDVCF | \
    grep -w $SNPID | \
    grep -v "^\#\#" | \
    cut -f 2,10- | \
    sed 's/:[0-9|\.]*\t/\t/g' | \
    sed 's/0|0/0/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | sed 's/1|1/2/g' | \
    datamash transpose | sort -k 1,1 > ${SCRATCHDIR}/snp_genotypes.tab

# Extract STR genotype and sample names to a file
bcftools query -r ${CHROM}:${POS}-${POS} -f '[%SAMPLE\t%TGT\n]' \
    ${IMPUTEDVCF} | sed 's/\//\t/' | sed 's/|/\t/' | \
    awk '{print $1 "\t" length($2) "\t" length($3)}' | \
    sort -k1,1 \
    > ${SCRATCHDIR}/str_genotypes.tab

# Extract phenotype, sex and sample names to a file
for bfile in $(cat ${COHORTFILES})
do
    cohort=$(basename $bfile | sed 's/-qc.bgs//')
    cat ${bfile}.fam | awk -v"cohort=$cohort" '{print $1"_"$2 "\t" $5 "\t" $6 "\t" cohort}'
done | sort -k1,1 > ${SCRATCHDIR}/phenotypes.tab

# Extract covariates and sample names to a file
cat ${COVARFILE} | awk -F " " '{print $1"_"$2 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $12 "\t" $18 "\t" $21}' | grep -v "FID IID" | \
    sort -k1,1 > ${SCRATCHDIR}/covars.tab

# Join all files
echo "sample,strgt1,strgt2,sex,phenotype,cohort,"${COVARCOLS}","${snpcols} | \
    sed 's/,/\t/g' > ${OUTDIR}/pgc_regression_data_med19.tab
join -t $'\t' ${SCRATCHDIR}/str_genotypes.tab ${SCRATCHDIR}/phenotypes.tab | \
    join -t $'\t' - ${SCRATCHDIR}/covars.tab | \
    join -t $'\t' - ${SCRATCHDIR}/snp_genotypes.tab \
    >> ${OUTDIR}/pgc_regression_data_med19.tab
#| \
#    join -t $'\t' - ${SCRATCHDIR}/snp_genotypes.tab \

