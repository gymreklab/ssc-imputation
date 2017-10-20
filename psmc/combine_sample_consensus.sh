#!/bin/bash

#SBATCH -A ddp268
#SBATCH -p shared
#SBATCH -t 1000
#SBATCH --get-user-env
#SBATCH --mem=8G

set -e -o pipefail

source params.sh

sample=$1

for chrom in $(seq 1 22)
do
    zcat ${chrompath}/${sample}.final.bam.fq_${chrom}.gz
done > ${samplepath}/${sample}.fq
gzip -f ${samplepath}/${sample}.fq

exit 0