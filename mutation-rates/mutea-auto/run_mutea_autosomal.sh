#!/bin/bash

source params.sh

CHROM=$1

if [ "x${SLURM_ARRAY_TASK_ID}" == "x" ]
then
    BATCH=$2
else
    BATCH=${SLURM_ARRAY_TASK_ID}
fi

DEBUGARGS=$3

# Get batch number padded with 0s to length 5
xBATCH=$(printf "%05d" ${BATCH})

python ${MUTEADIR}/mutea-auto/main_autosomal.py \
    --asdhet ${ASDT_VCF}/hipstr.chr${CHROM}.asdt.vcf.gz --vcf \
    --out ${OUTDIR}/batch_estimates/ssc_hipstr_mutea_chrom${CHROM}_batch${BATCH}.tab \
    --use-likelihoods \
    --loci ${OUTDIR}/batches/${CHROM}.${xBATCH} \
    --min_samples ${MINSAMPLES} \
    --min_mu ${MINMU} --max_mu ${MAXMU} \
    --min_beta ${MINBETA} --max_beta ${MAXBETA} \
    --min_pgeom ${MINPGEOM} --max_pgeom ${MAXPGEOM} \
    --stderrs fisher ${DEBUGARGS}