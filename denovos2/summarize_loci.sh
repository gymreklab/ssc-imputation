#!/bin/bash

source params.sh

#################### Update summary files ###########################
for chrom in $(seq 1 22)
do
    cmd="./summarize_loci.py \
	--allmutations ${FINALOUTDIR}/denovos_bylength.all_mutations_filtered.tab.gz \
	--loci ${FINALOUTDIR}/denovos_by_length_loci.bed \
	--annotations ${ANNDIR}/denovo_annotations.bed \
	--out ${OUTDIR}/tmp/denovos_bylength.locus_summary_${chrom}.bed \
	--output-mutations ${OUTDIR}/tmp/denovos_bylength.all_mutations_pass_${chrom}.tab \
	--filter-both-kids \
	--max-filtered-families ${MAXFILTFAM} \
	--min-children ${MINCHILDREN} \
	--chrom ${chrom} \
	--pthresh ${PTHRESH} > ${LOGDIR}/chr${chrom}.summ.out 2>&1"
    echo $cmd
done | xargs -P22 -I% -n1 sh -c "%"

cat ${OUTDIR}/tmp/denovos_bylength.locus_summary_1.bed | head -n 1 > ${FINALOUTDIR}/denovos_bylength.locus_summary.bed
cat ${OUTDIR}/tmp/denovos_bylength.all_mutations_pass_1.tab | head -n 1 > ${FINALOUTDIR}/denovos_bylength.all_mutations_pass.tab
for chrom in $(seq 1 22)
do
    cat ${OUTDIR}/tmp/denovos_bylength.locus_summary_${chrom}.bed | grep -v "^#" >> \
	${FINALOUTDIR}/denovos_bylength.locus_summary.bed
    cat ${OUTDIR}/tmp/denovos_bylength.all_mutations_pass_${chrom}.tab | grep -v "^#" >> \
	${FINALOUTDIR}/denovos_bylength.all_mutations_pass.tab
done

bgzip -f ${FINALOUTDIR}/denovos_bylength.locus_summary.bed
tabix -p bed -f ${FINALOUTDIR}/denovos_bylength.locus_summary.bed.gz
bgzip -f ${FINALOUTDIR}/denovos_bylength.all_mutations_pass.tab
tabix -b 2 -e -2 -f ${FINALOUTDIR}/denovos_bylength.all_mutations_pass.tab.gz

