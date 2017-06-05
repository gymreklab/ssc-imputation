#!/bin/bash

set -e

BAMFILE=$1
CHROM=$2
OUTBUCKET=$3
OUTFILE=$4
REFFA=$5

# Proceed
out=$(dirname ${OUTFILE})/$(basename ${OUTFILE} .gz)_${CHROM}.gz
sudo samtools mpileup -r ${CHROM} -C50 -uf $REFFA $BAMFILE | bcftools call -c - \
      | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > ${out}
aws s3 cp ${out} ${OUTBUCKET}/$(basename ${out})
rm ${out}
