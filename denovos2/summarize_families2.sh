#!/bin/bash

source params.sh

PERIOD=$1
AFTERFILTER=$2

if [ x$AFTERFILTER == x0 ]
then
    DIR=$OUTDIR
else
    DIR=$FINALOUTDIR
fi

tmpdir=$(mktemp -d)

echo $tmpdir

echo "family,child,status,numutations,prior_sum,posterior_sum,pthresh" | sed 's/,/\t/g' > ${DIR}/denovos_bylength_bychild2_period${PERIOD}.tab

echo "all mutation count $PERIOD"
pthresh=0
if [ x$AFTERFILTER == x0 ]
then
    cat ${DIR}/denovos_chr*_bylength.all_mutations.tab
else
    zcat ${DIR}/denovos_bylength.all_mutations_filtered.tab.gz
fi | grep -v chrom | \
    grep -v nan | \
    awk -v"period=$PERIOD" '($3==period)' | \
    sort -T ${OUTDIR}/tmp -k5,5 -k6,6 | \
    datamash -g 5,6,7 count 5 sum 4 sum 8 | awk -v"pthresh=${pthresh}" '{print $0 "\t" pthresh}' \
    >> ${DIR}/denovos_bylength_bychild2_period${PERIOD}.tab

for pthresh in 0.5 0.6 0.7 0.8 0.9 0.95 0.99 1
do
    echo ${pthresh} ${PERIOD}
    if [ x$AFTERFILTER == x0 ]
    then
	cat ${DIR}/denovos_chr*_bylength.all_mutations.tab
    else
	zcat ${DIR}/denovos_bylength.all_mutations_filtered.tab.gz
    fi | grep -v chrom | \
	grep -v nan | \
	awk -v"period=$PERIOD" -v"pthresh=${pthresh}" '($3==period && $8>=pthresh)' \
	> ${tmpdir}/ssc${pthresh}.tab
    cat ${tmpdir}/ssc${pthresh}.tab | \
	sort -T ${OUTDIR}/tmp -k5,5 -k6,6 | \
	datamash -g 5,6,7 count 5 sum 4 sum 8 | awk -v"pthresh=${pthresh}" '{print $0 "\t" pthresh}' \
	>> ${DIR}/denovos_bylength_bychild2_period${PERIOD}.tab
done 

