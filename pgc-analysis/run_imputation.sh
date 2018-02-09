#!/bin/bash

source params.sh

# Script below runs a single cohort/chrom
#for chrom in $(seq 1 22)
#do
chrom=21
    for cohort in $(cat ${COHORTFILES} | head -n 3)
    do
	cat impute_cohort_template.sh | sed "s/\$2/${chrom}/" | sed "s~\$1~${cohort}~" > run.sh
	chmod +x run.sh
	qsub -e $cohort_$chrom.err -o $cohort_$chrom.out ./run.sh
    done
#done
