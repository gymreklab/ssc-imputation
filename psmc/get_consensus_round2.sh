#!/bin/bash

set -e -o pipefail

BAMFILE=$1
CHROM=$2
OUTBUCKET=$3
OUTFILE=$4
REFFA=$5
LAB_AWS_ACCESS_KEY=$6
LAB_AWS_SECRET_KEY=$7
SSC_AWS_ACCESS_KEY=$8
SSC_AWS_SECRET_KEY=$9

# Proceed
out=$(dirname ${OUTFILE})/$(basename ${OUTFILE} .gz)_${CHROM}.gz
# Temporarily set default profile to ssc
export AWS_DEFAULT_PROFILE=ssc
export AWS_ACCESS_KEY_ID=${SSC_AWS_ACCESS_KEY}
export AWS_SECRET_ACCESS_KEY=${SSC_AWS_SECRET_KEY}
sudo -E samtools mpileup -r ${CHROM} -C50 -uf $REFFA $BAMFILE | bcftools call -c - \
      | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > ${out}
# Set back profile to default
export AWS_DEFAULT_PROFILE=default
export AWS_ACCESS_KEY_ID=${LAB_AWS_ACCESS_KEY}
export AWS_SECRET_ACCESS_KEY=${LAB_AWS_SECRET_KEY}
aws s3 cp ${out} ${OUTBUCKET}/$(basename ${out})
rm ${out}
