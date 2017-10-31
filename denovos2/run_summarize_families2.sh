#!/bin/bash

source params.sh
for period in $(seq 1 6)
do
	echo "./summarize_families2.sh ${period}"
done | xargs -P6 -I% -n1 sh -c "%"

cat ${FINALOUTDIR}/denovos_bylength_bychild2_period* | \
    head -n 1 > ${FINALOUTDIR}/denovos_bylength_bychild2.tab
cat ${FINALOUTDIR}/denovos_bylength_bychild2_period* | \
    awk '{print $NF "\t" $0}' | sort | \
    grep -v family | datamash -g 1,2,3,4 sum 5 sum 6 sum 7 mean 8 | \
    cut -f 1 --complement >> ${FINALOUTDIR}/denovos_bylength_bychild2.tab
