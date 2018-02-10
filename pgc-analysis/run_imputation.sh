#!/bin/bash

source params.sh

# Script below runs a single cohort/chrom
#for chrom in $(seq 1 22)
#do
chrom=17
    for cohort in $(cat ${COHORTFILES})
    do
	cat impute_cohort_template.sh | sed "s/\$2/${chrom}/" | sed "s~\$1~${cohort}~" > run.sh
	chmod +x run.sh
	qsub -e log_$cohort_$chrom.err -o log_$cohort_$chrom.out ./run.sh
    done
#done
