#!/bin/bash

OUTDIR=/storage/mgymrek/ssc-imputation/loo

for chrom in $(seq 1 22)
do
    VCF=/storage/mgymrek/ssc-imputation/filtered_vcfs/hipstr.chr${chrom}.allfilters.vcf.gz
    nohup ./compute_expected_metrics.py --vcf ${VCF} --out ${OUTDIR}/tmp/${chrom}.random.tab &
done
