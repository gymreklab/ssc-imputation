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
for chrom in $(seq 1 22)
do
    ./merge_cohort_results.sh ${chrom}
done

# Merge all files for each chrom - SNPs
for chrom in $(seq 1 22)
do
    echo $chrom
    ./merge_cohort_results_snps.sh ${chrom}
done

# Perform regression analysis
for chrom in $(seq 1 22)
do
    ./run_regression_wrapper.sh ${chrom}
#    ./pgc_regression.sh ${chrom} 2> logs/regr_log_${chrom}.err
done

