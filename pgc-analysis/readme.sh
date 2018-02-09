#!/bin/bash

# Get list of reference SNPs
for chrom in $(seq 1 22)
do
    ./get_ref_snps.sh ${chrom}
done

# Extract sample name, phenotype, sex, and covars to a single file
./extract_sample_info.sh

# Run imputation on all chroms/cohorts
./run_imputation.sh

# Merge all files for each chrom - TODO
./merge_cohort_results ${chrom}

# Perform regression analysis - TODO
./pgc_regression.sh ${chrom}
