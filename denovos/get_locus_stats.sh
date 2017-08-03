#!/bin/bash

source params.sh
OUTDIR=${BASEOUTDIR}/locus_stats

# Get HWE for each locus
for chrom in $(seq ${startchrom} ${endchrom})
do
    VCFFILE=/storage/s1saini/hipstr_genomewide/chr${chrom}/hipstr_calls_${chrom}.vcf.gz
#    /home/mgymrek/workspace/mgymrek-utils/vcf_hwe.py \
#	--vcf ${VCFFILE} \
#	--samples ${PARENTS} > ${OUTDIR}/ssc_hwe_chr${chrom}.tab
done

# Get number of de novo calls
for chrom in $(seq ${startchrom} ${endchrom})
do
    DFFILE=${BASEOUTDIR}/denovofinder/denovofinder_chr${chrom}.vcf.gz
#    zcat ${DFFILE} | grep -v "^#" | cut -f 1,2,10- | sed 's/\./0/g' | \
#	sed -E -e 's/SSC[0-9]*,SSC[0-9]*:-[0-9]*:-[0-9]*:-[0-9]*,-[0-9]*:-[0-9]*,-[0-9]*\t/1\t/g' | \
#	awk '{for(i=3;i<=NF;i++) t+=$i; print $1 "\t" $2 "\t" t; t=0}' > ${OUTDIR}/ssc_denovostats_chr${chrom}.tab
done


