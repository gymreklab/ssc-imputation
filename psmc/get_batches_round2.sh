#!/bin/bash

BATCHDIR=$1
BATCHSIZE=$2
DONE=$3

# Clear existing batch dir
mkdir -p ${BATCHDIR}
rm ${BATCHDIR}/*

# Batch files should have format bampath:chrom

BAMFILES=../metadata/ssc_parent_bampaths.txt
#DONE=completed_round1.txt
JOBSFILE=jobs_round2.txt

for bamfile in $(cat ${BAMFILES})
do
    for chrom in $(seq 1 22)
    do
	bam=$(basename $bamfile)
	outfile=$bam.fq_${chrom}.gz
	indone=$(grep -w $outfile $DONE)
	if [[ -z "${indone// }" ]]
	then
	    echo ${bamfile}:${chrom}
	else
	    continue
	fi
    done
done > ${JOBSFILE}

split -l ${BATCHSIZE} ${JOBSFILE} ${BATCHDIR}/bamfiles2
aws s3 sync ${BATCHDIR}/ s3://ssc-psmc/${BATCHDIR}
