#!/bin/bash

source params.sh

infile=${OUTDIR}/genome-wide/ssc_mutea_auto_unfiltered.bed.gz
outfile=${OUTDIR}/genome-wide/ssc_mutea_auto_scaled.bed

./filter_and_scale.py \
    --infile ${infile} \
    --outfile ${outfile} \
    --scale ${SCALE} \
    --gamma ${GAMMA}
bgzip -f ${OUTDIR}/genome-wide/ssc_mutea_auto_scaled.bed
tabix -p bed ${OUTDIR}/genome-wide/ssc_mutea_auto_scaled.bed.gz