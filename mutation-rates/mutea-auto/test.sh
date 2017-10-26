#!/bin/bash

sbatch -A ddp268 \
    -t 2000 \
    --mem=5G \
    -p shared \
    --job-name=testmutea \
    --get-user-env \
    -o testmutea_%a.out \
    -e testmutea_%a.err \
    --array=18,26,200 \
    ./run_mutea_autosomal.sh 2
