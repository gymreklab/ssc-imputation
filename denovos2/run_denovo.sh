#!/bin/bash

source params.sh

for chrom in $(seq ${startchrom} ${endchrom})
do
    echo "/home/mgymrek/workspace/STRDenovoTools/src/STRDenovoTools \
	--strvcf /storage/s1saini/hipstr_genomewide/chr${chrom}/hipstr_calls_${chrom}.vcf.gz \
	--fam /home/mgymrek/workspace/ssc-imputation/denovos/pedigree.fam \
	--max-num-alleles ${MAXALLELES} \
	--require-all-children \
	--require-num-children 2 \
	--min-coverage ${MINCOV} \
	--min-score ${MINSCORE} \
	--min-span-coverage ${MINSPANCOV} \
	--min-supp-reads ${MINSUPPREADS} \
	--posterior-threshold ${PTHRESH} \
	--out ${OUTDIR}/denovos_chr${chrom}_byseq"
done | xargs -P22 -I% -n1 sh -c "%"
