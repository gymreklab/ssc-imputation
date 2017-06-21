#!/bin/bash

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
	./process_denovo.sh ${chrom} ${family} ${father} ${mother} ${affected} ${unaffected}
    done
    cat ${OUTDIR}/bychrom/SSC_denovo_${family}_*.tab.filtered | head -n 1 > ${OUTDIR}/byfamily/${family}_denovos.tab
    cat ${OUTDIR}/bychrom/SSC_denovo_${family}_*.tab.filtered | grep -v pos >> ${OUTDIR}/byfamily/${family}_denovos.tab

    # Summarize
    ./summarize_family.py \
	--infile ${OUTDIR}/byfamily/${family}_denovos.tab \
	--outfile ${OUTDIR}/byfamily/${family}_denovos_summary.tab \
	--prob ${PROBTHRESH} --family ${family}
done
