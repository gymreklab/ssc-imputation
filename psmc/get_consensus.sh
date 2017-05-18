#!/bin/bash

set -e

BAMFILE=$1
OUTBUCKET=$2
OUTFILE=$3
REFFA=$4

# First, check if output file exists on s3
x=$(aws s3 ls ${OUTBUCKET}/$(basename ${OUTFILE}) | awk '{print $NF}')
test -z $x || exit 0

# Proceed
sudo samtools mpileup -C50 -uf $REFFA $BAMFILE | bcftools call -c - \
      | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > ${OUTFILE}
aws s3 cp ${OUTFILE} ${OUTBUCKET}/$(basename ${OUTFILE})
rm ${OUTFILE}
