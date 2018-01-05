#!/bin/bash

source params.sh

if [ "x${SLURM_ARRAY_TASK_ID}" == "x" ]
then
    CHROM=$1
else
    CHROM=${SLURM_ARRAY_TASK_ID}
fi

vcf=${FINALVCFS}/hipstr.chr${CHROM}.allfilters.vcf.gz

bcftools query -f "[%SAMPLE\t%GT\n]" ${vcf} | grep -v "\." | cut -f 1 | \
    sort | uniq -c > ${SAMPDIR}/sample_stats_${CHROM}.tab