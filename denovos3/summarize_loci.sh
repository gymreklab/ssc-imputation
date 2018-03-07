#!/bin/bash
#SBATCH -A csd568
#SBATCH --get-user-env 
#SBATCH --job-name Summarize_SSC_Denovo_Loci
#SBATCH -p shared
#SBATCH -t 1000
#SBATCH --mem=16gb
#SBATCH -o /oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered_mutations/logs/Summarize_SSC_Denovo_Loci.log
#SBATCH -e /oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered_mutations/logs/Summarize_SSC_Denovo_Loci.log         #SBATCH --mail-user=ilmitra@ucsd.edu
#SBATCH --mail-type=ALL
#SBATCH --array=1-21

date;hostname;pwd

source /home/ileena/ssc-imputation/denovos3/params.sh

#################### Update summary files ###########################
chrom=${SLURM_ARRAY_TASK_ID}
/home/ileena/ssc-imputation/denovos3/summarize_loci.py \
	--allmutations ${FINALOUTDIR}/denovos_bylength.all_mutations_filtered.tab.gz \
	--loci ${FINALOUTDIR}/denovos_by_length_loci.bed \
	--annotations ${ANNDIR}/denovo_annotations.bed \
	--out ${OUTDIR}/denovos_chr${chrom}_bylength.locus_summary.tab \
	--output-mutations ${OUTDIR}/denovos_chr${chrom}_bylength.all_mutations.tab \
	--filter-both-kids \
	--max-filtered-families ${MAXFILTFAM} \
	--min-children ${MINCHILDREN} \
	--chrom ${chrom} \
	--pthresh ${PTHRESH} > ${LOGDIR}/chr${chrom}.summ.out 2>&1





