#!/bin/bash

source params.sh

for chrom in $(seq ${startchrom} ${endchrom})
do
    echo "/home/mgymrek/workspace/STRDenovoTools/src/STRDenovoTools \
	--strvcf /storage/s1saini/hipstr_rerun/chr${chrom}/hipstr.chr${chrom}.with.1kg.vcf.gz \
	--fam /home/mgymrek/workspace/ssc-imputation/denovos/pedigree.fam \
	--max-num-alleles ${MAXALLELES} \
	--require-all-children \
	--require-num-children 2 \
	--min-coverage ${MINCOV} \
	--min-score ${MINSCORE} \
	--min-span-coverage ${MINSPANCOV} \
	--min-supp-reads ${MINSUPPREADS} \
	--posterior-threshold ${PTHRESH} \
	--out ${OUTDIR}/denovos_chr${chrom}_byseq > chr${chrom}_seq.out 2>&1"
done | xargs -P22 -I% -n1 sh -c "%"
