#!/bin/bash

CHROM=$1

cat run_imputation_cl.sh | sed "s~\$1~${CHROM}~" > run.sh
chmod +x run.sh
qsub -N impute${CHROM} -e logs/log_${CHROM}.err -o logs/log_${CHROM}.out ./run.sh
