#!/bin/bash

#PBS -lwalltime=2:00:00 -lnodes=1:ppn=1
#PBS -d /home/gymrek/workspace/ssc-imputation/pgc-analysis

source params.sh

chrom=$1

for cohort in $(cat ${COHORTFILES})
do
    ./get_pgc_vcf.sh ${cohort} ${chrom}
    ./impute_pgc.sh ${cohort} ${chrom}
done
