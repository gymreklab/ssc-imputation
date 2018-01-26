# Each example outputs columns:
# locus1       locus2  allele	freq	KL	r2	pval

# Get LD for a single SNP/STR pair
./snp_str_ld_calculator.py \
  --str-vcf /storage/s1saini/hipstr_rerun/chr5/hipstr.chr5.with.1kg.filtered.vcf.gz \
  --snp-vcf /storage/resources/datasets/SSC_SNP_v2/shapeit.chr5.with.ref.vcf.gz \
  --str-locus 5:153681115 \
  --snp-locus-rsid rs11740474 \
  --use-info-start \
  --max-dist 5000

./snp_str_ld_calculator.py \
  --str-vcf /storage/s1saini/hipstr_rerun/chr5/hipstr.chr5.with.1kg.filtered.vcf.gz \
  --snp-vcf /storage/resources/datasets/SSC_SNP_v2/shapeit.chr5.with.ref.vcf.gz \
  --str-locus 5:153681115 \
  --snp-locus 5:153680747 \
  --use-info-start \
  --max-dist 5000

# Get SNP-STR LD for all SNPs within some window of each STR
./snp_str_ld_calculator.py \
  --str-vcf /storage/s1saini/hipstr_rerun/chr21/hipstr.chr21.with.1kg.filtered.vcf.gz \
  --snp-vcf /storage/resources/datasets/SSC_SNP_v2/shapeit.chr21.with.ref.vcf.gz \
  --pairwise-snpstr \
  --max-dist 10000

# Get r2 for two different STR callsets (e.g. true vs. imputed)
./snp_str_ld_calculator.py \
  --str-vcf hipstr.1kg.EUR.filtered.vcf.gz \
  --str-vcf2 1kg.EUR.wgs.imputed.vcf.gz \
  --mincount 3

# Get per-allele r2 (can be used in all cases above as well)
./snp_str_ld_calculator.py \
  --str-vcf hipstr.1kg.EUR.filtered.vcf.gz \
  --str-vcf2 1kg.EUR.wgs.imputed.vcf.gz \
  --mincount 3 \
  --allele-r2

# Read pairs of STR/SNP to get LD for from a file (use either rsid or SNP locus)
# See help message for format of the loci file
# Also example of subsetting to set of certain samples (e.g. only parents)
# --user-info-start matches on INFO["START"] rather than on "POS" which HipSTR can change
./snp_str_ld_calculator.py \
  --str-vcf /storage/s1saini/hipstr_rerun/chr17/hipstr.chr17.with.1kg.filtered.vcf.gz \
  --snp-vcf /storage/resources/datasets/SSC_SNP_v2/shapeit.chr17.with.ref.vcf.gz \
  --loci-file-rsid /storage/mgymrek/ssc-imputation/pgc/bychrom/ldfile_17.tab \
  --samples ~/workspace/ssc-imputation/metadata/ssc_parent_ids.txt \
  --use-info-start --allele-r2
