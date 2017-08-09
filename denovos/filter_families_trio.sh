#!/bin/bash

set -e -o pipefail

source params.sh

# For each family
for family in $(cat ${FAMFILE} | cut -f 1 -d' '| sort | uniq)
do
    famfile=${BASEOUTDIR}/byfamily/SSC_denovo_${family}.tab
    echo ${famfile}
    ./filter_denovos_trio.py \
	--infile ${famfile} \
	--outfile ${famfile}.filtered \
	--minq-child ${MINQ} --minq-father ${MINQ} --minq-mother ${MINQ} \
	--minst-child ${MINST} --minst-father ${MINST} --minst-mother ${MINST} \
	--min-spanning-child ${MINSPAN} --min-spanning-father ${MINSPAN} --min-spanning-mother ${MINSPAN} \
	--min-percsupp-child ${MINSUPP} --min-percsupp-father ${MINSUPP} --min-percsupp-mother ${MINSUPP} \
	--parent-allele-count ${PARCOUNT} \
	--unit
done
