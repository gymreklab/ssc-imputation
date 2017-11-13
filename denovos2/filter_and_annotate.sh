#!/bin/bash

source params.sh

tmpdir=$(mktemp -d)

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
# Annotations:

# Clean list of annotations
rm -f ${ANNDIR}/annotations/*.txt

# Get list of all loci
cat ${OUTDIR}/denovos_chr*_bylength.locus_summary.tab | \
    grep -v chrom | cut -f 1-3 | \
    sort -k1,1 -k2,2n | uniq > ${ANNDIR}/all_loci.bed

# Gene annotations
cat ${GENES} | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
    sort -k1,1 -k2,2n | uniq | \
    datamash -g 1,2,3 unique 4 > ${ANNDIR}/annotations_gene.txt

# ExAC
cat ${EXAC} | grep -v chr | awk '{print $3 "\t" $5 "\t" $6 "\t" $18}' | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' |  \
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
# UTRs
cat ${UTR3} | grep -v track | sed 's/^chr//' | cut -f 1-3 | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -c -wa > \
    ${ANNDIR}/annotations_3utr.txt
cat ${UTR5} | grep -v track | sed 's/^chr//' | cut -f 1-3 | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -c -wa > \
    ${ANNDIR}/annotations_5utr.txt
# Promoters
cat ${PROMOTER1KB} | grep -v track | sed 's/^chr//' | cut -f 1-3 | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -c -wa > \
    ${ANNDIR}/annotations_promoter1kb.txt
cat ${PROMOTER3KB} | grep -v track | sed 's/^chr//' | cut -f 1-3 | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -c -wa > \
    ${ANNDIR}/annotations_promoter3kb.txt
cat ${PROMOTER5KB} | grep -v track | sed 's/^chr//' | cut -f 1-3 | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -c -wa > \
    ${ANNDIR}/annotations_promoter5kb.txt

# RNA binding sites 
cat ${RNABP} | sed 's/^chr//' | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
    sort -k1,1 -k2,2n | uniq | \
    datamash -g 1,2,3 sum 4 > ${ANNDIR}/annotations_rnabp.txt

# TF binding sites - TODO

# Known autism genes
for genescore in S 1 2 3 4 5 6
do
    cat ${ASDGENES} | awk -F',' -v"genescore=$genescore" '($(NF-2)==genescore)' | \
	cut -d',' -f 2 | sort > ${tmpdir}/${genescore}_sfarigenes.txt
    cat ${GENES} | grep -w -f ${tmpdir}/${genescore}_sfarigenes.txt - | grep -v "-" | \
	intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
	awk -v"genescore=$genescore" '{print $1 "\t" $2 "\t" $3 "\t" (($NF==".")?0:genescore)}' | \
	sort -k1,1 -k2,2n | uniq | \
	datamash -g 1,2,3 unique 4 > ${ANNDIR}/annotations_sfari${genescore}.txt
done
# Known autism CNVs-TODO

# ESTRs-TODO

# STR constraint scores
cat ${CONSTRAINT} | grep -v chrom | awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
    sort -k1,1 -k2,2n | uniq | \
    datamash -g 1,2,3 unique 4 > ${ANNDIR}/annotations_strconsZ.txt    
cat ${CONSTRAINT} | grep -v chrom | awk '{print $1 "\t" $2 "\t" $3 "\t" $(NF-1)}' | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
    sort -k1,1 -k2,2n | uniq | \
    datamash -g 1,2,3 unique 4 > ${ANNDIR}/annotations_strconsDIFF.txt    

# Motif
zcat ${HIPPROP} | grep -v chrom | cut -f 1-4 | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
    sort -k1,1 -k2,2n | uniq | \
    datamash -g 1,2,3 unique 4 > ${ANNDIR}/annotations_motif.txt

# Combine annotations
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
