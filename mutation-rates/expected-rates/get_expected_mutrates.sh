#!/bin/bash

source params.sh

tmpdir=$(mktemp -d)

for period in $(seq 1 6)
do
    logmu=$(cat ${MUTPARAMS} | awk -v"period=$period" '($1==period)' | awk '{print $2}')
    slope=$(cat ${MUTPARAMS} | awk -v"period=$period" '($1==period)' | awk '{print $3}')
    meanlen=$(cat ${MUTPARAMS} | awk -v"period=$period" '($1==period)' | awk '{print $4}')
    cat ${HIPREF} | awk -v"period=$period" '($4==period)' | \
	awk -v"logmu=$logmu" -v "slope=$slope" -v "meanlen=$meanlen" \
	'{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" ($3-$2+1-meanlen)*slope+logmu}' | \
	awk -v"maxval=$MAXVAL" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" ($5>maxval?maxval:$5)}'
done > ${tmpdir}/mutpred.bed

sort -k1,1 -k2,2n ${tmpdir}/mutpred.bed > ${OUTDIR}/predicted_str_mutrates_GRCh37.bed
