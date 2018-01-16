#!/bin/bash

#SBATCH -A csd568
#SBATCH --get-user-env
#SBATCH -t 5
#SBATCH -p shared

source params.sh

plink --vcf ${SNPS} \
    --chr ${CHROM} \
    --r2 \
    --maf ${MAF} \
    --out ${DIR}/tmp/tmp \
    --ld-window 1000000 \
    --ld-window-kb ${WINDOWKB}

echo "distance,r2" | sed 's/,/\t/g' > ${DIR}/snp_pairwise_ld.tab
cat ${DIR}/tmp/tmp.ld | grep -v CHR | awk '{print ($5-$2) "\t" $NF}' \
    >> ${DIR}/snp_pairwise_ld_chr${CHROM}.tab