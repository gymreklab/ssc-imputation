#!/bin/bash

source params.sh

# Separate BAM based on phase
./separate_bams.sh

# Run HipSTR separately on each phase in haploid mode
./run_hipstr_phased.sh

# Identify best SNP for each STR
./identify_snpstrs.sh

# Extract 10X + panel gts for each snpstr
./extract_phased_gts.sh
