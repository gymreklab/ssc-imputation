#!/bin/bash

source params.sh

# TODO: how to handle posterior=nan?

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

#################### Update summary files ###########################
# TODO remove calls if mutation in both kids
./summarize_loci.py \
    --allmutations ${FINALOUTDIR}/denovos_bylength.all_mutations_filtered.tab.gz \
    --loci ${FINALOUTDIR}/denovos_by_length_loci.bed \
    --annotations ${ANNDIR}/denovo_annotations.bed \
    --out ${FINALOUTDIR}/denovos_bylength.locus_summary.bed \
    --output-mutations ${FINALOUTDIR}/denovos_bylength.all_mutations_pass.tab \
    --filter-both-kids \
    --max-filtered-families ${MAXFILTFAM} \
    --pthresh ${PTHRESH}
bgzip -f ${FINALOUTDIR}/denovos_bylength.locus_summary.bed
tabix -p bed -f ${FINALOUTDIR}/denovos_bylength.locus_summary.bed.gz
