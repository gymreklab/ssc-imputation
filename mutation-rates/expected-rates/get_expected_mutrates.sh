#!/bin/bash

source params.sh

tmpdir=$(mktemp -d)
echo $tmpdir

for period in $(seq 1 6)
do
    logmu=$(cat ${MUTPARAMS} | awk -v"period=$period" '($1==period)' | awk '{print $2}')
    slope=$(cat ${MUTPARAMS} | awk -v"period=$period" '($1==period)' | awk '{print $3}')
    meanlen=$(cat ${MUTPARAMS} | awk -v"period=$period" '($1==period)' | awk '{print $4}')
    meanbeta=$(cat ${MUTPARAMS} | awk -v"period=$period" '($1==period)' | awk '{print $5}')
    meanp=$(cat ${MUTPARAMS} | awk -v"period=$period" '($1==period)' | awk '{print $6}')
    cat ${HIPREF} | awk -v"period=$period" '($4==period)' | \
	intersectBed -a stdin -b hipstr_unint_len.bed.gz -wa -wb | \
	awk -v"logmu=$logmu" -v "slope=$slope" -v "meanlen=$meanlen" -v "meanbeta=$meanbeta" -v "meanp=$meanp" \
	'{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" ($NF)*slope+logmu "\t" meanbeta "\t" meanp}' | \
	awk -v"maxval=$MAXVAL" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" ($5>maxval?maxval:$5) "\t" $6 "\t" $7}'
done > ${tmpdir}/mutpred.bed

sort -k1,1 -k2,2n ${tmpdir}/mutpred.bed | \
    intersectBed -a stdin -b ssc_central_alleles.bed.gz -sorted -wa -wb -loj | sed 's/\.$/0/' | \
    cut -f 1-7,11 > ${OUTDIR}/predicted_str_mutrates_GRCh37.bed
