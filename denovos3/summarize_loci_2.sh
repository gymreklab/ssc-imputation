#!/bin/bash
#SBATCH -A csd568
#SBATCH --get-user-env 
#SBATCH --job-name Summarize_SSC_Denovo_Loci
#SBATCH -p shared
#SBATCH -t 1000
#SBATCH --mem=16gb
#SBATCH -o /oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered_mutations/logs/Summarize_SSC_Denovo_Loci.log
#SBATCH -e /oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered_mutations/logs/Summarize_SSC_Denovo_Loci.log                                                          
#SBATCH --mail-user=ilmitra@ucsd.edu
#SBATCH --mail-type=ALL

date;hostname;pwd

source /home/ileena/ssc-imputation/denovos3/params.sh

#################### Update summary files ###########################

cat ${OUTDIR}/denovos_chr1_bylength.locus_summary.tab | head -n 1 > ${FINALOUTDIR}/denovos_bylength.locus_summary.bed
cat ${OUTDIR}/denovos_chr1_bylength.all_mutations.tab | head -n 1 > ${FINALOUTDIR}/denovos_bylength.all_mutations_pass.tab
for chrom in $(seq 1 22)
do
    cat ${OUTDIR}/denovos_chr${chrom}_bylength.locus_summary.tab | grep -v "^#" >> \
	${FINALOUTDIR}/denovos_bylength.locus_summary.bed
    cat ${OUTDIR}/denovos_chr${chrom}_bylength.all_mutations.tab | grep -v "^#" >> \
	${FINALOUTDIR}/denovos_bylength.all_mutations_pass.tab
done

bgzip -f ${FINALOUTDIR}/denovos_bylength.locus_summary.bed
tabix -p bed -f ${FINALOUTDIR}/denovos_bylength.locus_summary.bed.gz
bgzip -f ${FINALOUTDIR}/denovos_bylength.all_mutations_pass.tab
tabix -b 2 -e -2 -f ${FINALOUTDIR}/denovos_bylength.all_mutations_pass.tab.gz
