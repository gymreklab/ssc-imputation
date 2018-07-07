#!/bin/bash

source params.sh

HAPLOID="chr1"
for chrom in $(seq 2 22)
do
    HAPLOID="${HAPLOID},chr${chrom}"
done

for phase in 1 2
do
    /home/mgymrek/workspace/HipSTR/HipSTR \
	--bams ${OUTDIR}/NA12878_10x_phase${phase}.bam \
	--fasta ${REFFA} \
	--regions ${REGIONS} \
	--str-vcf ${OUTDIR}/NA12878_10x_phase${phase}.vcf.gz \
	--log ${OUTDIR}/NA12878_10x_phase${phase}.log.txt \
	--def-stutter-model --min-reads 10 \
	--use-unpaired \
	--haploid-chrs $HAPLOID
done

