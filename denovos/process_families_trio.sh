#!/bin/bash

set -e -o pipefail

source params.sh

# For each family
for family in $(cat ${FAMFILE} | cut -f 1 -d' '| sort | uniq)
do
    # Process all chroms
    father=$(cat ${FAMFILE} | awk -v"family=$family"  -F" " '(($1==family) && ($3!=0))' | cut -f 3 -d' ' | head -n 1)
    mother=$(cat ${FAMFILE} | awk -v"family=$family"  -F" " '(($1==family) && ($3!=0))' | cut -f 4 -d' ' | head -n 1)
    affected=$(cat ${FAMFILE} | awk -v"family=$family"  -F" " '(($1==family) && ($3!=0) && ($6==2))' | cut -f 2 -d' ' | head -n 1)
    unaffected=$(cat ${FAMFILE} | awk -v"family=$family"  -F" " '(($1==family) && ($3!=0) && ($6==1))' | cut -f 2 -d' ' | head -n 1)
    for chrom in $(seq $startchrom $endchrom)
    do
	echo ./process_denovo_trio.sh ${chrom} ${family} ${father} ${mother} ${affected} ${unaffected}
    done
done | xargs -P 4 -I% -n1 sh -c "%"

# Combine chroms
for family in $(cat ${FAMFILE} | cut -f 1 -d' '| sort | uniq)
do
    head -n 1 ${BASEOUTDIR}/bychrom/SSC_denovo_${family}_*.tab > ${BASEOUTDIR}/byfamily/SSC_denovo_${family}.tab
    cat ${BASEOUTDIR}/bychrom/SSC_denovo_${family}_*.tab | grep -v gtdata >> ${BASEOUTDIR}/byfamily/SSC_denovo_${family}.tab
    rm ${BASEOUTDIR}/bychrom/SSC_denovo_${family}_*.tab
done
