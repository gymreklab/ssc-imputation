#!/bin/bash

# Subset files to only have parents
./subset_vcfs.sh

# Get r2_imp per-locus
./get_r2_imp.sh

# Get r2_snp pairwise within window
./get_r2_snp.sh

# Phenotype simulations
./run_sims.sh # quant phenotype
./run_sims_cc.sh # case control
./run_examples.sh # examples for A, B plots
