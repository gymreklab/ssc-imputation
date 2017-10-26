#!/bin/bash

PERIOD=$1

source params.sh
tmpdir=$(mktemp -d)

echo $tmpdir

chrom="*"

echo "family,child,status,numutations,pthresh" | sed 's/,/\t/g' > ${OUTDIR}/denovos_bylength_bychild2_period${PERIOD}.tab

echo "all mutation count"
pthresh=0
cat ${OUTDIR}/denovos_chr${chrom}_bylength.all_mutations.tab | grep -v chrom | \
    awk -v"period=$PERIOD" '($3==period)' | \
    cut -f 4,5,6 | sort | uniq -c | \
    sed -e 's/ *//' -e 's/ /\t/' | \
    awk -v"pthresh=${pthresh}" '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" pthresh}'  \
    >> ${OUTDIR}/denovos_bylength_bychild2_period${PERIOD}.tab

for pthresh in 0.5 0.6 0.7 0.8 0.9 0.95 0.99 1
do
    echo ${pthresh}
    cat ${OUTDIR}/denovos_chr${chrom}_bylength.all_mutations.tab | grep -v chrom | \
	awk -v"period=$PERIOD" '($3==period)' | \
	awk -v"pthresh=${pthresh}" '($7>=pthresh)' > ${tmpdir}/ssc${pthresh}.tab
    cat ${tmpdir}/ssc${pthresh}.tab | \
	cut -f 4,5,6 | sort | uniq -c | \
	sed -e 's/ *//' -e 's/ /\t/' | \
	awk -v"pthresh=${pthresh}" '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" pthresh}'  \
	>> ${OUTDIR}/denovos_bylength_bychild2_period${PERIOD}.tab
done 

