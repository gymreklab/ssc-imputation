#!/bin/bash

#family=14679
#father=SSC12787
#mother=SSC12783
#proband=SSC12779
#sibling=SSC12870
family=12481
father=SSC04945
mother=SSC04943
proband=SSC04941
sibling=SSC04946

samtools view -H s3://sscwgs/$family/BAM/Sample_$father/analysis/$father.final.bam  > ${family}_${father}_father.sam
samtools view -H s3://sscwgs/$family/BAM/Sample_$mother/analysis/$mother.final.bam  > ${family}_${mother}_mother.sam
samtools view -H s3://sscwgs/$family/BAM/Sample_$proband/analysis/$proband.final.bam  > ${family}_${proband}_proband.sam
samtools view -H s3://sscwgs/$family/BAM/Sample_$sibling/analysis/$sibling.final.bam  > ${family}_${sibling}_sibling.sam

while IFS='' read -r line || [[ -n "$line" ]]; do
    chrom=$(echo $line | cut -f 1 -d' ')
    start=$(echo $line | cut -f 2 -d' ')
    end=$(echo $line | cut -f 3 -d' ')
    locus=$(echo $line | cut -f 4 -d' ')
    echo $locus
    samtools view s3://sscwgs/$family/BAM/Sample_$father/analysis/$father.final.bam ${chrom}:${start}-${end} >> ${family}_${father}_father.sam
    samtools view s3://sscwgs/$family/BAM/Sample_$mother/analysis/$mother.final.bam ${chrom}:${start}-${end} >> ${family}_${mother}_mother.sam
    samtools view s3://sscwgs/$family/BAM/Sample_$proband/analysis/$proband.final.bam ${chrom}:${start}-${end} >> ${family}_${proband}_proband.sam
    samtools view s3://sscwgs/$family/BAM/Sample_$sibling/analysis/$sibling.final.bam ${chrom}:${start}-${end} >> ${family}_${sibling}_sibling.sam
done < "loci.bed"

for f in $(ls ${family}*.sam)
do
    samtools view -bS $f | samtools sort - > $(basename $f .sam).bam
    samtools index $(basename $f .sam).bam
done

rm ${family}*.sam
