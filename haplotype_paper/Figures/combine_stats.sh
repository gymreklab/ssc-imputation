#!/bin/bash

TMPPATH=/storage/mgymrek/ssc-imputation/tmp

cat ${TMPPATH}/locus_stats_* > ${TMPPATH}/locus_stats.tab
cat ${TMPPATH}/sample_stats_* > ${TMPPATH}/sample_stats.tab
cat ${TMPPATH}/het_stats_* > ${TMPPATH}/het_stats.tab
cat ${TMPPATH}/mend_stats_* > ${TMPPATH}/mend_stats.tab
