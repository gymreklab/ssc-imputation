#!/bin/bash

#SBATCH -A ddp268
#SBATCH -p shared
#SBATCH -t 1000
#SBATCH --get-user-env

source params.sh

SAMPLE=$1

# Get psmcfa
/home/mgymrek/workspace/psmc/utils/fq2psmcfa -q20 ${samplepath}/${SAMPLE}.fq.gz > ${TMPDIR}/${SAMPLE}.psmcfa

# Run psmc
/home/mgymrek/workspace/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -d -o ${OUTDIR}/${SAMPLE}.psmc ${TMPDIR}/${SAMPLE}.psmcfa
