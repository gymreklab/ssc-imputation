# Get SNP pairwise LD
./snp_ld.sh

# Get SNP-STR pairwise LD - chr21
./snp_str_ld_calculator.py \
  --str-vcf ${STRS} \
  --snp-vcf ${SNPS} \
  --pairwise-snpstr
