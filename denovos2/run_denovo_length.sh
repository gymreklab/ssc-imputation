#!/bin/bash

startchrom=19
endchrom=22
OUTDIR=/storage/mgymrek/ssc-denovos/denovos2/denovocalls

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
	--posterior-threshold 0.9 \
	--combine-alleles-by-length \
	--out ${OUTDIR}/denovos_chr${chrom}_bylength
done
