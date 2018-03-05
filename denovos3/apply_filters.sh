#!/bin/bash
#SBATCH -A csd568
#SBATCH --get-user-env 
#SBATCH --job-name Apply_Filters_Denovos
#SBATCH -p shared
#SBATCH -t 1000
#SBATCH -o /oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered/logs/Apply_Filter_SSC_Phase1_Denovos.log
#SBATCH -e /oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered/logs/Apply_Filter_SSC_Phase1_Denovos.log 
#SBATCH --mail-user=ilmitra@ucsd.edu
#SBATCH --mail-type=ALL



source ~/ssc-imputation/denovos3/params.sh

#################### Apply filters to mutation files################
echo "#"$(head -n 1 ${OUTDIR}/denovos_chr1_bylength.all_mutations.tab) | sed 's/ /\t/g'  > ${FINALOUTDIR}/denovos_bylength.all_mutations_filtered.tab
for chrom in $(seq 1 22)
do
    cat ${OUTDIR}/denovos_chr${chrom}_bylength.all_mutations.tab | grep -v chrom | \
	awk '{print $1 "\t" $2 "\t" $2+1 "\t" $0}' | \
	intersectBed -a stdin -b ${FILTERDIR}/denovo_locus_filters.bed -v | \
	grep -v -w -f ${FILTERDIR}/denovo_family_filters.txt | \
	cut -f 1-3 --complement \
	>> ${FINALOUTDIR}/denovos_bylength.all_mutations_filtered.tab
done
bgzip -f ${FINALOUTDIR}/denovos_bylength.all_mutations_filtered.tab
tabix -b 2 -e -2 -c '#' -f ${FINALOUTDIR}/denovos_bylength.all_mutations_filtered.tab.gz

# Get list of loci
zcat ${FINALOUTDIR}/denovos_bylength.all_mutations_filtered.tab.gz | grep -v chrom | \
    cut -f 1,2 | uniq > ${FINALOUTDIR}/denovos_by_length_loci.bed

