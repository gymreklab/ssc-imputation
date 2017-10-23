#!/bin/bash

source params.sh

for chrom in $(seq 1 22)
do
    # How many batches?
    numbatches=$(ls ${OUTDIR}/batches/${chrom}.* | wc -l)
    sbatch -A ddp268 -t 2800 --array=0-"$(($numbatches-1))" \
	-p shared \
	--job-name=mutea_ssc_chr${chrom} --get-user-env \
	-o ${LOGDIR}/mutea_ssc_chr${chrom}.out -e ${LOGDIR}/mutea_ssc_chr${chrom}.err
	./run_mutea_autosomal.sh ${chrom}
done