#!/bin/bash

source params.sh

#for chrom in $(seq 1 22)
#do
chrom=20
    cat ${BESTSNPDIR}/chr${chrom}_snp_str_ld_50k.txt | \
	datamash -H groupby 1 max 7 -f | cut -f 1,2,7 > \
	${OUTDIR}/snpstrs/snpstrs_${chrom}.tab
#done
