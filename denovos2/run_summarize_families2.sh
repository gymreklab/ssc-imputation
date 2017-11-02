#!/bin/bash

source params.sh

AFTERFILTER=$1

if [ x$AFTERFILTER == x0 ]
then
    DIR=$OUTDIR
else
    DIR=$FINALOUTDIR
fi

for period in $(seq 1 6)
do
	echo "./summarize_families2.sh ${period}" $AFTERFILTER
done | xargs -P6 -I% -n1 sh -c "%"

cat ${DIR}/denovos_bylength_bychild2_period* | \
    head -n 1 > ${DIR}/denovos_bylength_bychild2.tab
cat ${DIR}/denovos_bylength_bychild2_period* | \
    awk '{print $NF "\t" $0}' | sort | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8}' | \
    grep -v family | datamash -g 1,2,3,4 sum 5 sum 6 sum 7 mean 8 | \
    cut -f 1 --complement >> ${DIR}/denovos_bylength_bychild2.tab
