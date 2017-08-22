#!/bin/bash

source params.sh

for chrom in $(seq ${startchrom} ${endchrom})
do
    /home/mgymrek/workspace/STRDenovoTools/src/STRDenovoTools \
	--strvcf /storage/s1saini/hipstr_genomewide/chr${chrom}/hipstr_calls_${chrom}.vcf.gz \
	--fam /home/mgymrek/workspace/ssc-imputation/denovos/pedigree.fam \
	--max-num-alleles 25 \
	--require-all-children \
	--require-num-children 2 \
	--min-coverage 10 \
	--min-score 0.9 \
	--posterior-threshold ${PTHRESH} \
	--out ${OUTDIR}/denovos_chr${chrom}
done
