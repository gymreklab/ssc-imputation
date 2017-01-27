#!/bin/bash

GATK=/storage/resources/source/GenomeAnalysisTK.jar
REFFA=/storage/resources/dbase/human/hs37d5/hs37d5.fa

# For testing
TESTGVCF1=/storage/s1saini/VCF/gvcf/11002/Sample_SSC03070/analysis/SSC03070.haplotypeCalls.raw.gvcf.vcf.gz
TESTGVCF2=/storage/s1saini/VCF/gvcf/11002/Sample_SSC03078/analysis/SSC03078.haplotypeCalls.er.raw.vcf.ava.gz
OUTDIR=/storage/mgymrek/ssc-imputation/gatk-hc

# http://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode

# CombineGVCFs - recommend batches of ~200 gVCFs into one

# Joint genotyping - GenotypeGVCFs
java -jar ${GATK} \
    -T GenotypeGVCFs \
    -R ${REFFA} \
    --variant ${TESTGVCF1} \
    --variant ${TESTGVCF2} \
    -o ${OUTDIR}/test.vcf

# Classic GATK Best Practices
