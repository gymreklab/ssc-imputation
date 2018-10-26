#!/bin/bash

source params.sh

# Get recombination hotspots
cat ${RECOMB} | awk '($5>=10)' | bedtools sort -i stdin > ${OUTDIR}/hotspots_hg19.bed

# For each STR, get distance from nearest hotspot
#cat ${HIPREF} | bedtools sort -i stdin | \
#    closestBed -a stdin -b ${OUTDIR}/hotspots_hg19.bed -d | \
#    cut -f 1,2,3,4,5,10,11 > ${OUTDIR}/str_dist_hotspots.bed

# For each STR, get local recombination rate
cat ${RECOMB} | awk '{print $1 "\t" $2-5000 "\t" $2+5000 "\t" $5}' | \
    intersectBed -a ${HIPREF} -b stdin -wa -wb | cut -f 1,2,3,4,5,9 > ${OUTDIR}/str_recomb_rate.bed

# In Jupyter:
# Density of STRs vs. dist from hotspot
# Mut rate vs. dist from hotspot
# Constraint vs. dist
# Heterozygosity
# Average length
