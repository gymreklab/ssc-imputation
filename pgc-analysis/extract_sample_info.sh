#!/bin/bash

source params.sh

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
echo "sample,sex,phenotype,cohort,"${COVARCOLS} | \
    sed 's/,/\t/g' > ${OUTDIR}/pgc_sample_data.tab
join -t $'\t' ${SCRATCHDIR}/phenotypes.tab ${SCRATCHDIR}/covars.tab \
    >> ${OUTDIR}/pgc_regression_data.tab

