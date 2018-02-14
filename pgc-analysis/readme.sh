#!/bin/bash

# Extract sample name, phenotype, sex, and covars to a single file
./extract_sample_info.sh

# Get files needed to make VCFs from plnk
for chrom in $(seq 1 22)
do
    ./get_ref_snps.sh ${chrom}
done

# Run imputation on all chroms/cohorts
for chrom in $(seq 1 22)
do
    ./run_imputation_cl_wrapper.sh $chrom
done

# Merge all files for each chrom
./merge_cohort_results.sh ${chrom}

# Perform regression analysis
./pgc_regression.sh ${chrom}
