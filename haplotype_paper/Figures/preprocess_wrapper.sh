#!/bin/bash

VCFPATH=/storage/s1saini/hipstr_rerun/
TMPPATH=/storage/mgymrek/ssc-imputation/tmp
FAMFILE=~/workspace/ssc-imputation/metadata/ssc_family_ids.txt

for chrom in $(seq 1 22)
do
    cmd="./preprocess_vcfs.sh ${VCFPATH} ${TMPPATH} ${FAMFILE} ${chrom}"
    echo $cmd
done | xargs -P 5 -I% -n 1 sh -c "%"
