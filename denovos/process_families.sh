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
	echo "[process_families.sh] Running process_denovo.sh..."
	./process_denovo.sh ${chrom} ${family} ${father} ${mother} ${affected} ${unaffected}
	echo "[process_families.sh] Running filter_denovos.py..."
	./filter_denovos.py \
	    --infile ${OUTDIR}/bychrom/SSC_denovo_${family}_${chrom}.tab \
	    --outfile ${OUTDIR}/bychrom/SSC_denovo_${family}_${chrom}.tab.filtered \
	    --minq-child ${MINQ} --minq-father ${MINQ} --minq-mother ${MINQ} \
	    --minst-child ${MINST} --minst-father ${MINST} --minst-mother ${MINST} \
	    --min-spanning-child ${MINSPAN} --min-spanning-father ${MINSPAN} --min-spanning-mother ${MINSPAN} \
	    --min-percsupp-child ${MINSUPP} --min-percsupp-father ${MINSUPP} --min-percsupp-mother ${MINSUPP} \
	    --parent-allele-count ${PARCOUNT} \
	    --unit
    done
    cat ${OUTDIR}/bychrom/SSC_denovo_${family}_*.tab.filtered | head -n 1 > ${OUTDIR}/byfamily/${family}_denovos.tab
    cat ${OUTDIR}/bychrom/SSC_denovo_${family}_*.tab.filtered | grep -v pos >> ${OUTDIR}/byfamily/${family}_denovos.tab

    # Summarize
    ./summarize_family.py \
	--infile ${OUTDIR}/byfamily/${family}_denovos.tab \
	--outfile ${OUTDIR}/byfamily/${family}_denovos_summary.tab \
	--prob ${PROBTHRESH} --family ${family}
done
