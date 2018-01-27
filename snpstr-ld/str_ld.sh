#!/bin/bash

#SBATCH -A csd568
#SBATCH --get-user-env
#SBATCH -t 2800
#SBATCH -p shared

source params.sh

./snp_str_ld_calculator.py \
    --str-vcf ${STRS} \
    --snp-vcf ${SNPS} \
    --max-dist ${MAXDIST} \
    --min-maf ${MAF} \
    --usefilter \
    --region 21:15449764-48119754 \
    --pairwise-snpstr #> ${DIR}/snp_str_ld_chr21.tab 
