#!/bin/bash

source params.sh

tmpdir=$(mktemp -d -p ${TMPLOC})

# Clean existing batch data
rm -f ${OUTDIR}/batches/*

for chrom in $(seq 1 22)
do
    echo "Getting batches for chrom ${chrom}..."
    cat ${HIPREF} | awk -v"chrom=$chrom" '($1==chrom)' > ${tmpdir}/loc_${chrom}.bed
    # Split chromosome into batches
    split -l ${BATCHSIZE} -d -a 5 ${tmpdir}/loc_${chrom}.bed ${OUTDIR}/batches/${chrom}"."
done