#!/bin/bash

#PBS -lwalltime=2:00:00 -lnodes=1:ppn=1

source params.sh

#COHORT=$1
#CHROM=$2
COHORT=/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/scz/wave2/v1/cobg_dir_genome_wide/scz_aber_eur-qc.bgs
CHROM=21

./get_pgc_vcf.sh ${COHORT} ${CHROM}
./impute_pgc.sh ${COHORT} ${CHROM}

# Remove intermediate files
rm ${OUTDIR}/imputed/chr${CHROM}/tmp/$(basename $COHORT)*