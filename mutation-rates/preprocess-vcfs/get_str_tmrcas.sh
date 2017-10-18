#!/bin/bash

#SBATCH -A ddp268
#SBATCH -t 1000
#SBATCH -p shared
#SBATCH --get-user-env
#SBATCH -e str_tmrcas.err
#SBATCH -o str_tmrcas.out

source params.sh

# Process each sample individually
#for rawpsmc in $(ls ${RAW_PSMC}/*.psmc)
#do
#    sample=$(basename ${rawpsmc} .psmc)
#    ./get_str_tmrcas_bysample.sh ${sample}
#done

# Merge by chromosome
for chrom in $(seq 1 22)
do
    cat ${TMPDIR}/*.strpsmc.bed | awk -v"chrom=$chrom" '($1==chrom)' | sort -k1,1 -k2,2n | bgzip -c > ${OUT_PSMC_STR}/str_tmrcas_chr${chrom}.bed.gz
    tabix -f -p bed ${OUT_PSMC_STR}/str_tmrcas_chr${chrom}.bed.gz
done
