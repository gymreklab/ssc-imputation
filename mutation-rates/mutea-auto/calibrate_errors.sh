#!/bin/bash

source params.sh

./calibrate_errors.py \
    --ests ${OUTDIR}/test/ssc_hipstr_mutea_codis.tab \
    --truthnp codis_truth_np.bed \
    --scale ${SCALE}