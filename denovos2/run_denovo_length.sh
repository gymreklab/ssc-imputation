#!/bin/bash

source params.sh

for chrom in $(seq ${startchrom} ${endchrom})
do
    echo "/home/mgymrek/workspace/STRDenovoTools/src/STRDenovoTools \
	--strvcf /storage/s1saini/hipstr_rerun/chr${chrom}/hipstr.chr${chrom}.with.1kg.filtered.vcf.gz \
	--fam /home/mgymrek/workspace/ssc-imputation/denovos/pedigree.fam \
	--max-num-alleles ${MAXALLELES} \
	--require-all-children \
	--require-num-children 2 \
	--min-coverage ${MINCOV} \
	--min-score ${MINSCORE} \
	--min-span-coverage ${MINSPANCOV} \
	--min-supp-reads ${MINSUPPREADS} \
	--posterior-threshold ${PTHRESH} \
	--combine-alleles-by-length \
        --output-all-loci \
        --mutation-models /storage/mgymrek/ssc-denovos/mutea-results/predicted_str_mutrates_GRCh37.bed \
	--out ${OUTDIR}/denovos_chr${chrom}_bylength > ${LOGDIR}/chr${chrom}.out 2>&1"
done | xargs -P22 -I% -n1 sh -c "%"
