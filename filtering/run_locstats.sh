#!/bin/bash

source params.sh

sbatch \
    -A csd568 \
    --array=1-22 \
    -p shared \
    --mem=3G \
    -t 2000 \
    --get-user-env \
    --job-name=filter_hipstr_vcfs \
    -o ${LOGDIR}/locstats_%a.out -e ${LOGDIR}/locstats_%a.err \
    ./get_locstats.sh