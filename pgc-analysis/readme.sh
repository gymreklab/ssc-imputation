#!/bin/bash

# TODO make wrappers to run all chroms/cohorts

# Get list of reference SNPs
for chrom in $(seq 1 22)
do
    ./get_ref_snps.sh ${chrom}
done

# Extract sample name, phenotype, sex, and covars to a single file
./extract_sample_info.sh

# Script below runs a single cohort/chrom
# TODO make wrapper to submit for all cohorts/chroms. sub $1/$2
# see http://geneticcluster.org/pages/7/Tutorial.html
./impute_cohort_template.sh
