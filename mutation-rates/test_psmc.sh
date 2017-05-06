#!/bin/bash

# Following instructions here https://github.com/lh3/psmc

BAM=/storage/resources/datasets/gtex/bams/SRR2167793.bam
REFFA=/storage/resources/dbase/human/hg19/Homo_sapiens_assembly19.fasta

# Get whole genome consensus sequence
#samtools mpileup -C50 -uf $REFFA $BAM | bcftools call -c - \
#      | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > diploid.fq.gz

# Get psmcfa
#/home/mgymrek/workspace/psmc/utils/fq2psmcfa -q20 diploid.fq.gz > diploid.psmcfa

# Run psmc
#/home/mgymrek/workspace/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -d -o diploid.psmc diploid.psmcfa
