#!/bin/bash

source params.sh

tmpdir=$(mktemp -d -p ${TMPLOC})

for chrom in $(seq 1 22)
do
    echo "Getting batches for chrom ${chrom}..."
    cat ${HIPREF} | grep "^${chrom}" > ${tmpdir}/loc_${chrom}.bed
    # Split chromosome into batches
    split -l ${BATCHSIZE} -d -a 5 ${tmpdir}/loc_${chrom}.bed ${OUTDIR}/batches/${chrom}"."
    # Split stutter file into batches
    echo "Getting stutter batches for chrom ${chrom}..."
    for f in $(ls ${OUTDIR}/batches/${chrom}.*)
    do
	batchfile=${OUTDIR}/batches/$(basename $f)
	batch=$(basename $f | cut -d'.' -f 2)
	echo ${batchfile} ${STUTTERFILE}
	intersectBed -a ${STUTTERFILE} -b ${batchfile} -wa -u -f 1 > ${OUTDIR}/stutterfiles/ssc_hipstr_stutter_chr${chrom}_batch${batch}.bed
    done
done