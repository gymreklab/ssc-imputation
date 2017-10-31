#!/bin/bash

source params.sh

#################### Update summary files ###########################
./summarize_loci.py \
    --allmutations ${FINALOUTDIR}/denovos_bylength.all_mutations_filtered.tab.gz \
    --loci ${FINALOUTDIR}/denovos_by_length_loci.bed \
    --annotations ${ANNDIR}/denovo_annotations.bed \
    --out ${FINALOUTDIR}/denovos_bylength.locus_summary.bed \
    --output-mutations ${FINALOUTDIR}/denovos_bylength.all_mutations_pass.tab \
    --filter-both-kids \
    --max-filtered-families ${MAXFILTFAM} \
    --min-children ${MINCHILDREN} \
    --pthresh ${PTHRESH}
bgzip -f ${FINALOUTDIR}/denovos_bylength.locus_summary.bed
tabix -p bed -f ${FINALOUTDIR}/denovos_bylength.locus_summary.bed.gz
bgzip -f ${FINALOUTDIR}/denovos_bylength.all_mutations_pass.tab
tabix -b 2 -e -2 -f ${FINALOUTDIR}/denovos_bylength.all_mutations_pass.tab.gz

