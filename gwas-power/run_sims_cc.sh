#!/bin/bash

source params.sh

~/workspace/str-association-tests/power_analysis/snpstr_power_simulator.py \
    --snp-vcf ${OUTDIR}/phasedsnps_chr21_parents.vcf.gz \
    --str-vcf ${OUTDIR}/hipstr_chr21_parents.vcf.gz \
    --imputed-str-vcf ${OUTDIR}/hipstr_chr21_parents_imputed.vcf.gz \
    --beta 0.1 \
    --case-control \
    --numsim 100 \
    --use-best-snp ${OUTDIR}/r2_snp.tab \
    --outprefix ${OUTDIR}/snpstr_power_sims_casecontrol.tab

