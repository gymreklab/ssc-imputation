#!/bin/bash

source params.sh

SAMPLE=$1

raw=${RAW_PSMC}/${SAMPLE}.psmc
out=${TMPDIR}/${SAMPLE}.strpsmc.bed

echo "Processing PSMC for sample ${SAMPLE} ${raw} ${out}"

# Get weighted average of all intervals overlapping each STR locus
# NOTE!! looks like PSMC numbers are missing last two digits? binned by 100?

# mean genome-wide het
meanhet=$(cat ${raw} | awk '($1=="DC")' | awk '{print ($4-$3) "\t" ($4-$3)*$6}' | datamash sum 1 sum 2 | awk '{print $2/$1}')
scale=$(echo "${GWAVG}/${meanhet}" | bc -l)
echo "${sample} ${meanhet} ${scale}" > ${LOGDIR}/${SAMPLE}.strpsmc.log

cat ${raw} | awk '($1=="DC")' | \
    awk '{print $2 "\t" $3"00" "\t" $4"99" "\t" $6}' | \
    intersectBed -a ${HIPREF} -b stdin -wb -wa | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" ($9-$8) "\t" ($9-$8)*$10}' | \
    sort -k 1,1 -k2,2n | \
    datamash -g1,2,3 sum 4 sum 5 | \
    awk -v"sample=$SAMPLE" -v"scale=${scale}" '{print $1 "\t" $2 "\t" $3 "\t" ($5/$4)*scale "\t" sample}' > ${out}

