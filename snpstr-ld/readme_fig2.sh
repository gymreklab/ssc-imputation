#!/bin/bash

# Allele-r2 for SSC LOO
./snp_str_ld_calculator.py \
    --str-vcf /storage/s1saini/hipstr_allfilters/phased_feb18/hipstr.chr21.phased.vcf.gz \
    --str-vcf2 /storage/s1saini/manuscript_strsnp/fig3/loo/l1o.imputed.str.vcf.gz \
    --mincount 3 --usefilter \
    --samples /home/mgymrek/workspace/ssc-imputation/metadata/ssc_parent_ids.txt \
    --allele-r2 > /storage/mgymrek/ssc-imputation/loo/l1o_alleler2.tab

# Allele-r2 for SSC onekg
./snp_str_ld_calculator.py \
    --str-vcf /storage/s1saini/manuscript_strsnp/fig3/1kg/EUR/hipstr.1kg.EUR.filtered.vcf.gz \
    --str-vcf2 /storage/s1saini/manuscript_strsnp/fig3/1kg/EUR/1kg.EUR.wgs.imputed.vcf.gz \ 
    --mincount 3 --usefilter \
    --allele-r2 > /storage/mgymrek/ssc-imputation/onekg/onekg_alleler2.tab
