#!/bin/bash

source params.sh

# TODO get these from a file
slope1=0.05
slope2=0.056
slope3=0.056
slope4=0.067
slope5=0.06
slope6=0.06

FORCEBETA=0.2

echo "period,logmu,slope,meanlen,meanbeta,meanp" | sed 's/,/\t/g' > ../expected-rates/rate_model_params.tab
# Get mean vals by period
for period in $(seq 1 6)
do
    if [ x$period == x1 ]
    then
	slope=$slope1
    elif [ x$period == x2 ]
    then
	slope=$slope2
    elif [ x$period == x3 ]
    then
	slope=$slope3
    elif [ x$period == x4 ]
    then
	slope=$slope4
    elif [ x$period == x5 ]
    then
	slope=$slope5
    elif [ x$period == x6 ]
    then
	slope=$slope6
    else
	slope="NA"
    fi
    cat ${HIPREF} | awk -v"period=$period" '($4==period)' | \
	intersectBed -a ${OUTDIR}/genome-wide/ssc_mutea_auto_scaled.bed.gz -b stdin | \
	awk '($7>0)' | \
	awk '{print ($3-$2+1) "\t" $0}' | \
	datamash mean 1 mean 5 mean 6 mean 7 | \
	awk -v "beta=$FORCEBETA" -v"period=$period" -v"slope=$slope" '{print period "\t" $2 "\t" slope "\t" $1 "\t" beta "\t" $4}'
done >> ../expected-rates/rate_model_params.tab

# Get central alleles per locus
zcat ${OUTDIR}/genome-wide/ssc_mutea_auto_unfiltered.bed.gz | cut -f 1-3,9 | bgzip -c > \
    ../expected-rates/ssc_central_alleles.bed.gz
