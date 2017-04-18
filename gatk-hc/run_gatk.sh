#!/bin/bash

GATK=/storage/resources/source/GenomeAnalysisTK.jar
REFFA=/storage/resources/dbase/human/hs37d5/hs37d5.fa
BUNDLE=/storage/resources/dbase/human/hs37d5/gatk_resource_bundle/

# For testing
TESTGVCF1=/storage/s1saini/VCF/gvcf/11002/Sample_SSC03070/analysis/SSC03070.haplotypeCalls.raw.gvcf.vcf.gz
TESTGVCF2=/storage/s1saini/VCF/gvcf/11002/Sample_SSC03078/analysis/SSC03078.haplotypeCalls.er.raw.vcf.ava.gz
OUTDIR=/storage/mgymrek/ssc-imputation/gatk-hc

# http://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode

# CombineGVCFs - recommend batches of ~200 gVCFs into one

# Joint genotyping - GenotypeGVCFs
#java -jar ${GATK} \
#    -T GenotypeGVCFs \
#    -R ${REFFA} \
#    --variant ${TESTGVCF1} \
#    --variant ${TESTGVCF2} \
#    -o ${OUTDIR}/test.vcf

# Classic GATK Best Practices
java -jar ${GATK} \
    -T VariantRecalibrator \
    -R ${REFFA} \
    -input ${OUTDIR}/test.vcf \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${BUNDLE}/hapmap_3.3.b37.vcf \
    -resource:omni,known=false,training=true,truth=false,prior=12.0 ${BUNDLE}/1000G_omni2.5.b37.vcf \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${BUNDLE}/1000G_phase1.snps.high_confidence.b37.vcf \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${BUNDLE}/dbsnp_138.b37.vcf \
    -mode SNP \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
    -recalFile ${OUTDIR}/test.recal \
    -tranchesFile ${OUTDIR}/test.tranche \
    -rscriptFile ${OUTDIR}/test.plots.R
