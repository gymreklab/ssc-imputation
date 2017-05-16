#!/bin/bash

BAMFILE=$1
OUTFILE=$2
OUTBUCKET=$3
REFFA=$4

samtools mpileup -C50 -uf $REFFA $BAMFILE | bcftools call -c - \
      | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > ${OUTFILE}
aws s3 cp ${OUTFILE} ${OUTBUCKET}/${OUTFILE}
