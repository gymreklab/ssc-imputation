#!/bin/bash

source params.sh

CHROM=$1

cat pgc_regression.sh | sed "s~\$1~${CHROM}~" > run.sh
chmod +x run.sh
qsub -e logs/log_regr_${CHROM}.err -o logs/log_regr_${CHROM}.out ./run.sh
