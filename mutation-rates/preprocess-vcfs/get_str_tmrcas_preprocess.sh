#!/bin/bash

#SBATCH -A ddp268
#SBATCH -t 2800
#SBATCH -p shared
#SBATCH --get-user-env
#SBATCH -e str_tmrcas.err
#SBATCH -o str_tmrcas.out

source params.sh

# Process each sample individually
for rawpsmc in $(ls ${RAW_PSMC}/*.psmc)
do
    sample=$(basename ${rawpsmc} .psmc)
    ./get_str_tmrcas_bysample.sh ${sample}
done
