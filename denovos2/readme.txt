# Run denovo calling
./run_denovo_length.sh

# Get filters and annotations
./filter_and_annotate.sh
./apply_filters.sh

# Summarize families
./run_summarize_families2.sh 1

# Summarize all loci. Use filters learned from previous step
./summarize_loci.sh
