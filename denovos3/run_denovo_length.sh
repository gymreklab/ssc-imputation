#!/bin/bash
#SBATCH -A csd568
#SBATCH --get-user-env 
#SBATCH --job-name SSC_Phase1_Denovos
#SBATCH -p shared
#SBATCH -t 1000
#SBATCH --array=1-21
#SBATCH -o /oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/logs/SSC_Phase1_Denovos_%a.out
#SBATCH -e /oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/logs/SSC_Phase1_Denovos_%a.err



if [ "x${SLURM_ARRAY_TASK_ID}" == "x" ]
then
    chrom=$1
else
    chrom=${SLURM_ARRAY_TASK_ID}
fi

source /home/ileena/ssc-imputation/denovos3/params.sh

STRDenovoTools \
    --strvcf ${VCFDIR}/hipstr.chr${chrom}.allfilters.vcf.gz \
    --fam /home/ileena/ssc-imputation/denovos3/pedigree.fam \
    --max-num-alleles ${MAXALLELES} \
    --require-all-children \
    --require-num-children 2 \
    --min-coverage ${MINCOV} \
    --min-score ${MINSCORE} \
    --min-span-coverage ${MINSPANCOV} \
    --min-supp-reads ${MINSUPPREADS} \
    --posterior-threshold ${PTHRESH} \
    --combine-alleles-by-length --round-alleles --include-invariant \
    --output-all-loci \
    --mutation-model ${DATADIR}/predicted_str_mutrates_GRCh37.bed \
    --out ${OUTDIR}/denovos_chr${chrom}_bylength
