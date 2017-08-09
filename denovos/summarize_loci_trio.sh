#!/bin/bash

set -e -o pipefail

source params.sh

# Each call:
# - decide if affected/unaffected have length mutations
# - sum for each locus

cat ${BASEOUTDIR}/byfamily/*.filtered | grep -v chrom | \
    awk -F"\t" -v"probthresh=$PROBTHRESH" \
    '{print $1 "\t" $2 "\t" $4 "\t" \
         (($19=="False" && $11=="False" && $6>probthresh)?1:0) "\t" \
         (($20=="False" && $12=="False" && $9>probthresh)?1:0)}'| \
    sort -k1,1 -k2,2n | \
    datamash -g1,2,3 sum 4 count 4 sum 5 count 5 > ${BASEOUTDIR}/per_locus_summary.txt
