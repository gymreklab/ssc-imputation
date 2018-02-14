#!/bin/bash
source params.sh

EXLOCUS=21:23934706
echo "case control"
~/workspace/str-association-tests/power_analysis/snpstr_power_simulator.py \
    --snp-vcf ${OUTDIR}/phasedsnps_chr21_parents.vcf.gz \
    --str-vcf ${OUTDIR}/hipstr_chr21_parents.vcf.gz \
    --imputed-str-vcf ${OUTDIR}/hipstr_chr21_parents_imputed.vcf.gz \
    --beta 0.7 \
    --numsim 1 \
    --case-control \
    --str-region ${EXLOCUS} \
    --write-phenotypes --write-genotypes \
    --outprefix ${OUTDIR}/snpstr_power_sims_cc_example

# Fig. 4c
echo "quant"
~/workspace/str-association-tests/power_analysis/snpstr_power_simulator.py \
    --snp-vcf ${OUTDIR}/phasedsnps_chr21_parents.vcf.gz \
    --str-vcf ${OUTDIR}/hipstr_chr21_parents.vcf.gz \
    --imputed-str-vcf ${OUTDIR}/hipstr_chr21_parents_imputed.vcf.gz \
    --beta 0.3 \
    --numsim 1 \
    --str-region ${EXLOCUS} \
    --write-phenotypes --write-genotypes \
    --outprefix ${OUTDIR}/snpstr_power_sims_quant_example

