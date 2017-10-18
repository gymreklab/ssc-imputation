#!/bin/bash

source params.sh

for sample in $(ls -l $chrompath | grep mgymrek | awk '{print $NF}' | cut -f 1 -d'.' | uniq)
do
    sbatch \
	-J ${sample} \
	-o ${OUTDIR}/log/${sample}.con.out \
	-e ${OUTDIR}/log/${sample}.con.err \
	./combine_sample_consensus.sh ${sample}
    sbatch \
	-J ${sample} \
	-d singleton \
	-o ${OUTDIR}/log/${sample}.psmc.out \
	-e ${OUTDIR}/log/${sample}.psmc.err \
	./psmc_from_consensus.sh ${sample}
done