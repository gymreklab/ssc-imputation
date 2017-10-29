#!/bin/bash

source params.sh

for chrom in $(seq 1 22)
do
    echo $chrom $(ls ${OUTDIR}/batches/${chrom}.* | wc -l) \
	$(cat ${OUTDIR}/batches/${chrom}.* | wc -l)
done