#!/bin/bash

source params.sh

# TODO get these from a file
int1=-5.15
slope1=0.051

int2=-5.17
slope2=0.047

int3=-5.62
slope3=0.051

int4=-5.58
slope4=0.046

int5=-5.42
slope5=0.035

int6=-5.24
slope6=0.025

FORCEBETA=0.01

echo "period,logmu,slope,meanlen,meanbeta,meanp" | sed 's/,/\t/g' > ../expected-rates/rate_model_params.tab
# Get mean vals by period
for period in $(seq 1 6)
do
    if [ x$period == x1 ]
    then
	int=$int1
	slope=$slope1
    elif [ x$period == x2 ]
    then
	int=$int2
	slope=$slope2
    elif [ x$period == x3 ]
    then
	int=$int3
	slope=$slope3
    elif [ x$period == x4 ]
    then
	int=$int4
	slope=$slope4
    elif [ x$period == x5 ]
    then
	int=$int5
	slope=$slope5
    elif [ x$period == x6 ]
    then
	int=$int6
	slope=$slope6
    else
	int="NA"
	slope="NA"
    fi
    cat ${HIPREF} | awk -v"period=$period" '($4==period)' | \
	intersectBed -a ${OUTDIR}/genome-wide/ssc_mutea_auto_scaled.bed.gz -b stdin | \
	awk '($7>0)' | \
	awk '{print ($3-$2+1) "\t" $0}' | \
	datamash mean 1 mean 5 mean 6 mean 7 | \
	awk -v "beta=$FORCEBETA" -v"period=$period" -v"slope=$slope" -v"intercept=$int" '{print period "\t" intercept "\t" slope "\t" "0" "\t" beta "\t" $4}'
done >> ../expected-rates/rate_model_params.tab

# Get central alleles per locus
zcat ${OUTDIR}/genome-wide/ssc_mutea_auto_unfiltered.bed.gz | cut -f 1-3,9 | bgzip -c > \
    ../expected-rates/ssc_central_alleles.bed.gz
# Get uninterrupted length per locus
zcat ${CONSTRAINTDIR}/GRCh37.hipstr_reference_sorted_properties.tab.gz | \
    grep -v chrom | cut -f 1-3,6 | \
    bgzip -c > \
    ../expected-rates/hipstr_unint_len.bed.gz