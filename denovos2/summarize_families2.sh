#!/bin/bash

PERIOD=$1

source params.sh
tmpdir=$(mktemp -d)

echo $tmpdir

echo "family,child,status,numutations,prior_sum,posterior_sum,pthresh" | sed 's/,/\t/g' > ${FINALOUTDIR}/denovos_bylength_bychild2_period${PERIOD}.tab

echo "all mutation count"
pthresh=0
zcat ${FINALOUTDIR}/denovos_bylength.all_mutations_filtered.tab.gz | grep -v chrom | \
    grep -v nan | \
    awk -v"period=$PERIOD" '($3==period)' | \
    sort -k5,5 -k6,6 | \
    datamash -g 5,6,7 count 5 sum 4 sum 8 | awk -v"pthresh=${pthresh}" '{print $0 "\t" pthresh}' \
    >> ${FINALOUTDIR}/denovos_bylength_bychild2_period${PERIOD}.tab

for pthresh in 0.5 0.6 0.7 0.8 0.9 0.95 0.99 1
do
    echo ${pthresh}
    zcat ${FINALOUTDIR}/denovos_bylength.all_mutations_filtered.tab.gz | grep -v chrom | \
	grep -v nan | \
	awk -v"period=$PERIOD" -v"pthresh=${pthresh}" '($3==period && $8>=pthresh)' \
	> ${tmpdir}/ssc${pthresh}.tab
    cat ${tmpdir}/ssc${pthresh}.tab | \
	sort -k5,5 -k6,6 | \
	datamash -g 5,6,7 count 5 sum 4 sum 8 | awk -v"pthresh=${pthresh}" '{print $0 "\t" pthresh}' \
	>> ${FINALOUTDIR}/denovos_bylength_bychild2_period${PERIOD}.tab
done 

