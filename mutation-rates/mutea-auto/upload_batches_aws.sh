#!/bin/bash

#SBATCH -A ddp268
#SBATCH -p shared
#SBATCH -t 10
#SBATCH --get-user-env

source params.sh

BATCHDIR=${OUTDIR}/batches_aws/
SUPERBATCHDIR=${OUTDIR}/superbatches_aws/

tmpdir=$(mktemp -d -p ${TMPLOC})

# Clean existing batch data
rm -f ${BATCHDIR}/*
rm -f ${SUPERBATCHDIR}/*

# Make individual batches
for chrom in $(seq 1 22)
do
    echo "Getting batches for chrom ${chrom}..."
    cat ${HIPREF} | awk -v"chrom=$chrom" '($1==chrom)' > ${tmpdir}/loc_${chrom}.bed
    # Split chromosome into batches
    split -l ${AWSBATCHSIZE} -d -a 5 ${tmpdir}/loc_${chrom}.bed ${BATCHDIR}/${chrom}"."
done

# Make super batches
ls -l ${BATCHDIR}/* | grep -v total | awk -F'/' '{print $NF}' > ${tmpdir}/batches.txt
split -l ${AWSSUPERBATCHSIZE} -d -a 5 ${tmpdir}/batches.txt ${SUPERBATCHDIR}/superbatch"."

aws s3 sync ${BATCHDIR} s3://ssc-mutea/batches/
aws s3 sync ${SUPERBATCHDIR} s3://ssc-mutea/superbatches/