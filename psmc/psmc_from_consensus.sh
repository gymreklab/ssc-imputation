#!/bin/bash

SAMPLE=$1
CONSENSUSDIR=/oasis/projects/nsf/ddp268/mgymrek/ssc-quads/psmc/consensus
REFFA=/oasis/projects/nsf/ddp268/mgymrek/dbase/Homo_sapiens_assembly19.fasta
OUTDIR=/oasis/projects/nsf/ddp268/mgymrek/ssc-quads/psmc/raw_psmc/
TMPDIR=/oasis/scratch/comet/$USER/temp_project/

# Get psmcfa
/home/mgymrek/workspace/psmc/utils/fq2psmcfa -q20 ${CONSENSUSDIR}/${SAMPLE}.fq.gz > ${TMPDIR}/${SAMPLE}.psmcfa

# Run psmc
/home/mgymrek/workspace/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -d -o ${OUTDIR}/${SAMPLE}.psmc ${TMPDIR}/${SAMPLE}.psmcfa
