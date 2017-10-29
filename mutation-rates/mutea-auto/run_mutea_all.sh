#!/bin/bash

source params.sh

for chrom in $(seq 1 22)
do
    # How many batches?
    numbatches=$(ls ${OUTDIR}/batches/${chrom}.* | wc -l)
    sbatch -A ddp268 --array=0-"$(($numbatches-1))" \
	-p shared \
	--mem=5G -t 1000 \
	--job-name=mutea_ssc_chr${chrom} --get-user-env \
	-o ${LOGDIR}/mutea_ssc_chr${chrom}_batch%a.out -e ${LOGDIR}/mutea_ssc_chr${chrom}_batch%a.err \
	./run_mutea_autosomal.sh ${chrom}
done