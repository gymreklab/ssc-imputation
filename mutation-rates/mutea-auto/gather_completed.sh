#!/bin/bash

source params.sh

cat ${OUTDIR}/batch_estimates/*.tab | \
    sort -k1,1 -k2,2n | uniq | bgzip -c > ${OUTDIR}/genome-wide/ssc_mutea_auto_unfiltered.bed.gz
tabix -f -p bed ${OUTDIR}/genome-wide/ssc_mutea_auto_unfiltered.bed.gz