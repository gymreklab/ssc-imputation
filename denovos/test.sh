#!/bin/bash

source params.sh

VCFFILE=/storage/mgymrek/ssc-denovos/test/hipstr_calls_12_145.vcf.gz
OUTDIR=/storage/mgymrek/ssc-denovos/test/
SNPVCF=/storage/mgymrek/ssc-denovos/test/shapeit.chr12.vcf.gz 

${DENOVOFINDER} \
    --str-vcf ${VCFFILE} \
    --fam ${FAMFILE} \
    --snp-vcf ${SNPVCF} \
    --denovo-vcf ${OUTDIR}/$(basename ${VCFFILE} .vcf.gz )_denovos_withsnps.vcf.gz

${DENOVOFINDER} \
    --str-vcf ${VCFFILE} \
    --fam ${FAMFILE} \
    --denovo-vcf ${OUTDIR}/$(basename ${VCFFILE} .vcf.gz )_denovos_trio.vcf.gz
