#!/bin/bash

source params.sh

for chrom in $(seq 1 22)
do
    aws s3 cp ${ASDT_VCF}/hipstr.chr${chrom}.asdt.vcf.gz s3://ssc-strvcf/
    aws s3 cp ${ASDT_VCF}/hipstr.chr${chrom}.asdt.vcf.gz.tbi s3://ssc-strvcf/
done
