#!/bin/bash

sample=$1
chrompath=/oasis/projects/nsf/ddp268/mgymrek/ssc-quads/psmc//consensus_by_chrom/

for chrom in $(seq 1 22)
do
    zcat ${chrompath}/${sample}.final.bam.fq_${chrom}.gz
done > ${sample}.fq
gzip ${sample}.fq