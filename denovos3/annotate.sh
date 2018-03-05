#!/bin/bash
#SBATCH -A csd568
#SBATCH --get-user-env 
#SBATCH --job-name Annotate_SSC_Phase1_Denovos
#SBATCH -p shared
#SBATCH -t 1000
#SBATCH -o /oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered_mutations/logs/Annotate_SSC_Phase1_Denovos.log
#SBATCH -e /oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered_mutations/logs/Annotate_SSC_Phase1_Denovos.log 
#SBATCH --mail-user=ilmitra@ucsd.edu
#SBATCH --mail-type=ALL

source ~/ssc-imputation/denovos3/params.sh

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
# Introns
cat ${INTRON} | grep -v track | sed 's/^chr//' | cut -f 1-3 | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -c -wa > \
    ${ANNDIR}/annotations_introns.txt

# Splice sites
cat ${DONOR} | sed 's/^chr//' | sort -k1,1 -k2,2n | \
    closestBed -a ${ANNDIR}/all_loci.bed -b stdin -D b | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
    sort -k1,1 -k2,2n | uniq > ${ANNDIR}/annotations_donor.txt
cat ${ACCEPTOR} | sed 's/^chr//' | sort -k1,1 -k2,2n | \
    closestBed -a ${ANNDIR}/all_loci.bed -b stdin -D b | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
    sort -k1,1 -k2,2n | uniq > ${ANNDIR}/annotations_acceptor.txt

# TF binding sites (GM12878)
cat ${TF} | sed 's/^chr//' | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
    sort -k1,1 -k2,2n | uniq | \
    datamash -g 1,2,3 sum 4 > ${ANNDIR}/annotations_tfbs.txt
for colnum in $(seq 6 126)
do
    tf=$(head -n 1 ${TFALL} | cut -f ${colnum} | sed 's/_GM12878//')
    cat ${TFALL} | cut -f 1-3,${colnum} | grep -v chrom | sed 's/^chr//' | \
	intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
	awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
	sort -k1,1 -k2,2n | uniq | \
	datamash -g 1,2,3 sum 4 > ${ANNDIR}/annotations_${tf}.txt
done

# Histone mods (GM12878)
for colnum in $(seq 6 14)
do
    hmod=$(head -n 1 ${HISTONE} | cut -f ${colnum} | sed 's/_GM12878//')
    cat ${HISTONE} | cut -f 1-3,${colnum} | grep -v chrom | sed 's/^chr//' | \
	intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
	awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
	sort -k1,1 -k2,2n | uniq | \
	datamash -g 1,2,3 sum 4 > ${ANNDIR}/annotations_${hmod}.txt
done

# RNA binding sites - all
cat ${RNABP} | sed 's/^chr//' | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
    sort -k1,1 -k2,2n | uniq | \
    datamash -g 1,2,3 sum 4 > ${ANNDIR}/annotations_rnabp.txt

# RNABP - specifically TARDBP
cat ${RNABP2} | sed 's/^chr//' | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $NF}' | \
    sort -k1,1 -k2,2n | uniq | \
    datamash -g 1,2,3 sum 4 > ${ANNDIR}/annotations_rnabpTARDBP.txt

# Known autism genes (within 10kb)
for genescore in S 1 2 3 4 5 6
do
    cat ${ASDGENES} | awk -F',' -v"genescore=$genescore" '($(NF-2)==genescore)' | \
	cut -d',' -f 2 | sort > ${tmpdir}/${genescore}_sfarigenes.txt
    cat ${GENES} | grep -w -f ${tmpdir}/${genescore}_sfarigenes.txt - | grep -v "-" | \
	awk '{print $1 "\t" $2-10000 "\t" $3+10000}' | \
	intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -c > \
	${ANNDIR}/annotations_sfari${genescore}.txt
done
# Known autism CNVs
cat ${ASDCNV} | grep -v chrom | awk '($3>$2)' | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -wa -wb -loj | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" ($NF=="."?"0":$NF)}' |  \
    sort -k1,1 -k2,2n | uniq | \
    datamash -g 1,2,3 sum 4 > ${ANNDIR}/annotations_sfaricnv.txt

# ESTRs
cat ${ESTRS} | awk '($6<0.05)' | awk '{print $2 "\t" $3 "\t" $3+1}' | sed 's/^chr//' | sed 's/\.0//' | \
    grep -v best | grep -v NA | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -c -wa > \
    ${ANNDIR}/annotations_estrs.txt
    
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

# Brain expressed tissues
cat ${BRAIN} | sed 's/chr//' | \
    intersectBed -a ${ANNDIR}/all_loci.bed -b stdin -c -wa > \
    ${ANNDIR}/annotations_brain.txt

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
