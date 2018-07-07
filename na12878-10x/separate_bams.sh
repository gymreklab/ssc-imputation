#!/bin/bash

source params.sh

for phase in 1 2
do
    samtools view -h $BAM | grep "^@\|HP:i:${phase}" | \
	samtools view -bS > ${OUTDIR}/NA12878_10x_phase${phase}.bam
    samtools index ${OUTDIR}/NA12878_10x_phase${phase}.bam
done

