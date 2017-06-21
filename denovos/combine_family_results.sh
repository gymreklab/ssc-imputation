#!/bin/bash

source params.sh

cat ${OUTDIR}/byfamily/*_denovos_summary.tab | \
    grep -v family | cut -f 1 --complement | \
    sort -k 1,1 -k 2,2n | \
    datamash -g 1,2 count 1 sum 3 sum 4 sum 5 sum 6 | \
    awk '{print $0 "\t" ($4+$5)/$3 "\t" ($6+$7)/$3}' | sort -k8 -r > \
    ${OUTDIR}/ssc_denovo_mutation_counts.tab
