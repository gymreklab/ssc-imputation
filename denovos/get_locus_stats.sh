#!/bin/bash

source params.sh
OUTDIR=${BASEOUTDIR}/locus_stats

# Get HWE for each locus
for chrom in $(seq ${startchrom} ${endchrom})
do
    VCFFILE=/storage/s1saini/hipstr_genomewide/chr${chrom}/hipstr_calls_${chrom}.vcf.gz
    /home/mgymrek/workspace/mgymrek-utils/vcf_hwe.py \
	--vcf ${VCFFILE} \
	--samples ${PARENTS} \
	--sim 1000 \
	> ${OUTDIR}/ssc_hwe_chr${chrom}.tab
done

# Get number of de novo calls
for chrom in $(seq ${startchrom} ${endchrom})
do
    DFFILE=${BASEOUTDIR}/denovofinder/denovofinder_chr${chrom}.vcf.gz
#    zcat ${DFFILE} | grep -v "^#" | cut -f 1,2,10- | sed 's/\./0/g' | \
#	sed -E -e 's/SSC[0-9]*,SSC[0-9]*:-[0-9]*:-[0-9]*:-[0-9]*,-[0-9]*:-[0-9]*,-[0-9]*\t/1\t/g' | \
#	awk '{for(i=3;i<=NF;i++) t+=$i; print $1 "\t" $2 "\t" t; t=0}' > ${OUTDIR}/ssc_denovostats_chr${chrom}.tab
done

# Get stutter params, call rate for each locus
for chrom in $(seq ${startchrom} ${endchrom})
do
    VCFFILE=/storage/s1saini/hipstr_genomewide/chr${chrom}/hipstr_calls_${chrom}.vcf.gz
#    echo "chrom,pos,inframe_pgeom,inframe_up,inframe_down,outframe_pgeom,outframe_up,outframe_down,an,end,dp,dsnp,nfilt" | sed 's/,/\t/g' > ${OUTDIR}/ssc_gtstats_chr${chrom}.tab
#    bcftools query -f '%CHROM\t%POS\t%INFRAME_PGEOM\t%INFRAME_UP\t%INFRAME_DOWN\t%OUTFRAME_PGEOM\t%OUTFRAME_UP\t%OUTFRAME_DOWN\t%AN\t%END\t%DP\t%DSNP\t%NFILT\n' ${VCFFILE} >> ${OUTDIR}/ssc_gtstats_chr${chrom}.tab
done
