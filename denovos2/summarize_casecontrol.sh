#!/bin/bash

source params.sh

echo "#### By length ####"

# Overall, get number of mutations in case/control
echo "all" $(cat ${OUTDIR}/denovos_chr*_bylength.locus_summary.tab | \
    grep -v chrom | \
    datamash sum 12 sum 13 sum 15 sum 16)

# By period
for period in $(seq 1 6)
do
    echo ${period} $(cat ${OUTDIR}/denovos_chr*_bylength.locus_summary.tab | \
	awk -v"period=${period}" '($4==period)' | \
	grep -v chrom | \
	datamash sum 12 sum 13 sum 15 sum 16)
done

# Coding
echo "coding" $(cat ${OUTDIR}/denovos_chr*_bylength.locus_summary.tab | grep -v chrom | \
    awk '{print "chr"$0}' | \
    intersectBed -a stdin -b ${CODING} -u | \
    datamash sum 12 sum 13 sum 15 sum 16)
    

echo "#### By seq ####"

# Overall, get number of mutations in case/control
echo "all" $(cat ${OUTDIR}/denovos_chr*_byseq.locus_summary.tab | \
    grep -v chrom | \
    datamash sum 12 sum 13 sum 15 sum 16)

# By period
for period in $(seq 1 6)
do
    echo ${period} $(cat ${OUTDIR}/denovos_chr*_byseq.locus_summary.tab | \
	awk -v"period=${period}" '($4==period)' | \
	grep -v chrom | \
	datamash sum 12 sum 13 sum 15 sum 16)
done

# Coding
echo "coding" $(cat ${OUTDIR}/denovos_chr*_byseq.locus_summary.tab | grep -v chrom | \
    awk '{print "chr"$0}' | \
    intersectBed -a stdin -b ${CODING} -u | \
    datamash sum 12 sum 13 sum 15 sum 16)
    
