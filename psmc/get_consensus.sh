#!/bin/bash

set -e

BAMFILE=$1
OUTBUCKET=$2
OUTFILE=$3
REFFA=$4

# First, check if output file exists on s3
#x=$(aws s3 ls ${OUTBUCKET}/$(basename ${OUTFILE}) | awk '{print $NF}')
#test -z $x || exit 0

# Proceed
for chrom in $(seq 1 22)
do
out=$(dirname ${OUTFILE})/$(basename ${OUTFILE} .gz)_${chrom}.gz
sudo samtools mpileup -r ${chrom} -C50 -uf $REFFA $BAMFILE | bcftools call -c - \
      | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > ${out}
aws s3 cp ${out} ${OUTBUCKET}/$(basename ${out})
rm ${out}
done
