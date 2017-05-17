#!/bin/bash

set -e

BAMFILE=$1
OUTBUCKET=$2
OUTFILE=$3
REFFA=$4

sudo samtools mpileup -C50 -uf $REFFA $BAMFILE | head -n 10000 | bcftools call -c - \
      | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > ${OUTFILE}
aws s3 cp ${OUTFILE} ${OUTBUCKET}/$(basename ${OUTFILE})
rm ${OUTFILE}
