#!/bin/bash

source params.sh

########################################
# Which loci to filter

# Clean loci filters
rm -f ${FILTERDIR}/locus_filter*.bed

# Segdup
cat ${OUTDIR}/denovos_chr*_bylength.locus_summary.tab | \
    grep -v chrom | cut -f 1-3 | awk '{print "chr"$0}' | \
    intersectBed -a stdin -b ${SEGDUP} -wa | sed 's/^chr//' | \
    sort -k1,1 -k2,2n | uniq > \
    ${FILTERDIR}/locus_filter_segdup.bed

# Too few families called
cat ${OUTDIR}/denovos_chr*_bylength.locus_summary.tab | \
    grep -v chrom | \
    awk -v"minchildren=$MINCHILDREN" -v"minchildrencol=$MINCHILDRENCOL" \
    '($minchildrencol<minchildren)' | cut -f 1-3 | \
    sort -k1,1 -k2,2n | uniq > \
    ${FILTERDIR}/locus_filter_minfamilies.bed

# Combine all loci filters
cat ${FILTERDIR}/locus_filter*.bed | sort -k1,1 -k2,2n | uniq > \
    ${FILTERDIR}/denovo_locus_filters.bed

########################################
# Which families to filter

# Clean family filters
rm -f ${FILTERDIR}/family_filter*.txt

# This family had way too many denovos called in the affected child
# Use sample IDs instead of family IDs, which are just numbers
echo "SSC02229" > ${FILTERDIR}/family_filters_outlier.txt
echo "SSC02241" >> ${FILTERDIR}/family_filters_outlier.txt

# Combine all family filters
cat ${FILTERDIR}/*.txt | sort | uniq > ${FILTERDIR}/denovo_family_filters.txt

########################################
# Which calls to remove - TODO

########################################
# Annotations:

# Clean list of annotations
rm -f ${ANNDIR}/annotations/*.txt

# Get list of all loci
cat ${OUTDIR}/denovos_chr*_bylength.locus_summary.tab | \
    grep -v chrom | cut -f 1-3 | \
    sort -k1,1 -k2,2n | uniq > ${ANNDIR}/all_loci.bed

# ExAC
cat ${EXAC} | grep -v chr | awk '{print $3 "\t" $5 "\t" $6 "\t" $2}' | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
    sort -k1,1 -k2,2n | uniq | \
    datamash -g 1,2,3 unique 4 > ${ANNDIR}/annotations_gene.txt
cat ${EXAC} | grep -v chr | awk '{print $3 "\t" $5 "\t" $6 "\t" $18}' | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
    sort -k1,1 -k2,2n | uniq | \
    datamash -g 1,2,3 unique 4 > ${ANNDIR}/annotations_ExACmisZ.txt
cat ${EXAC} | grep -v chr | awk '{print $3 "\t" $5 "\t" $6 "\t" $20}' | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
    sort -k1,1 -k2,2n | uniq | \
    datamash -g 1,2,3 unique 4 > ${ANNDIR}/annotations_ExACpLI.txt

# Coding
cat ${CODING} | grep -v track | sed 's/^chr//' | cut -f 1-3 | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -c -wa > \
    ${ANNDIR}/annotations_coding.txt

# STR constraint scores - TODO

# Combine annotations
tmpdir=$(mktemp -d)
cp ${ANNDIR}/all_loci.bed ${tmpdir}/all_loci.bed
header="chrom,start,end"
for f in $(ls ${ANNDIR}/annotations*.txt)
do
    annname=$(basename $f .txt | sed 's/annotations_//')
    header=${header}","${annname}
    cat $f | cut -f 1-3 --complement | \
	paste ${tmpdir}/all_loci.bed - > ${tmpdir}/del.bed
    mv ${tmpdir}/del.bed ${tmpdir}/all_loci.bed
done
echo ${header} | sed 's/,/\t/g' > ${ANNDIR}/denovo_annotations.bed
cat ${tmpdir}/all_loci.bed >> ${ANNDIR}/denovo_annotations.bed
