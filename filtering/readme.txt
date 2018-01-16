# Run HipSTR call level filters
./run_filter_calls.sh

# Get locus stats
./run_locstats.sh

# Combine locus stats and set VCF filters
./run_set_locus_filters.sh

# Get sample stats
./run_sampstats.sh

# Get Mendelian inheritance for chr21
sbatch ./vcf_mend.sh 21
