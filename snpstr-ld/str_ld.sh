#!/bin/bash

#SBATCH -A csd568
#SBATCH --get-user-env
#SBATCH -t 2800
#SBATCH -p shared

source params.sh

# Get STR positions
#zcat ${STRS} | grep -v "^#" | grep PASS | cut -f 2 > ${DIR}/tmp/str_positions.txt

# Get samples
bcftools query -l ${STRS} | grep -w -f ${PARENTS} > ${SAMPLES}
nsamp=$(cat ${SAMPLES} | wc -l)

# Loop through positions, get all SNPs, get LD with each
while IFS='' read -r line || [[ -n "$line" ]]; do
    strpos=$line
    # Extract SNPs
    start=$((strpos-${WINDOWKB}000))
    end=$((strpos+${WINDOWKB}000))
    bcftools query -S ${SAMPLES} -r ${CHROM}:${strpos}-${strpos} -f '[%GB\n]' ${STRS} | \
	sed 's/\./100000|100000/' | sed 's/|/\t/' | awk '{print $1 + $2}' | sed 's/200000/./' \
	> ${DIR}/tmp/strgt.tab
    bcftools query -S ${SAMPLES} -r ${CHROM}:${start}-${end} -f '%INFO/AC\t%POS[\t%GT]\n' ${SNPS} | \
	awk -v"minac=$MINAC" '($1>=minac)' | cut -f 1 --complement | \
	shuf | head -n ${MAXSNPS} | \
	sed 's/0|0/0/g' | sed 's/1|1/2/g' | sed 's/0|1/1/g' | sed 's/1|0/1/g' | \
	datamash transpose > ${DIR}/tmp/snpgt.tab
    nsnp=$(cat ${DIR}/tmp/snpgt.tab | awk '{print NF}' | head -n 1)
    for snpnum in $(seq 1 ${nsnp})
    do
	cat ${DIR}/tmp/snpgt.tab | cut -f ${snpnum} | awk '(NR>1)' | paste - ${DIR}/tmp/strgt.tab | \
	    grep -v "\." | datamash ppearson 1:2 mean 1 | \
	    awk -v"strpos=$strpos" -v"snppos=$(head -n 1 ${DIR}/tmp/snpgt.tab | cut -f ${snpnum})" -v"chrom=$CHROM" \
	    '{print chrom "\t" strpos "\t" snppos "\t" (snppos-strpos) "\t" $1*$1 "\t" $2/2}'
    done
done < ${DIR}/tmp/str_positions.txt > ${DIR}/snp_str_ld_chr${CHROM}.tab
