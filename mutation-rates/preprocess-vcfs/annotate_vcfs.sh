#!/bin/bash

source params.sh

if [ "x${SLURM_ARRAY_TASK_ID}" == "x" ]
then
    CHROM=$1
else
    CHROM=${SLURM_ARRAY_TASK_ID}
fi

psmcfile=${OUT_PSMC_STR}/str_tmrcas_chr${CHROM}.bed.gz
strvcf=${HIPSTR_CALLS}/hipstr.chr${CHROM}.with.1kg.filtered.vcf.gz
outvcf=${OUT_VCF}/hipstr.chr${CHROM}.asdt.vcf.gz

./annotate_vcf.py --invcf ${strvcf} --asdt ${psmcfile} | \
    bcftools annotate -S ${SAMPLES} - | \
    bgzip -c > ${outvcf}
tabix -f -p vcf ${outvcf}