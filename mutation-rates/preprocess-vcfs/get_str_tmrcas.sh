#!/bin/bash

source params.sh

if [ "x${SLURM_ARRAY_TASK_ID}" == "x" ]
then
    CHROM=$1
else
    CHROM=${SLURM_ARRAY_TASK_ID}
fi


# Merge by chromosome
cat ${TMPDIR}/*.strpsmc.bed | awk -v"chrom=$chrom" '($1==chrom)' | sort -k1,1 -k2,2n | bgzip -c > ${OUT_PSMC_STR}/str_tmrcas_chr${CHROM}.bed.gz
tabix -f -p bed ${OUT_PSMC_STR}/str_tmrcas_chr${CHROM}.bed.gz

