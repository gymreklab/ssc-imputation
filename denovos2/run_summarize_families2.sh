#!/bin/bash

source params.sh
for period in $(seq 1 6)
do
	echo "./summarize_families2.sh ${period}"
done | xargs -P6 -I% -n1 sh -c "%"

cat ${OUTDIR}/denovos_bylength_bychild2_period* | \
    head -n 1 > ${OUTDIR}/denovos_bylength_bychild2.tab
cat ${OUTDIR}/denovos_bylength_bychild2_period* | \
    awk '{print $5 "\t" $1 "\t" $2 "\t" $3 "\t" $4}' | sort | \
    grep -v family | datamash -g 1,2,3,4 sum 5 | \
    awk '{print $0 "\t" $1}' | cut -f 1 --complement >> ${OUTDIR}/denovos_bylength_bychild2.tab
