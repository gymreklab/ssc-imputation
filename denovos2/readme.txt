# Run denovo calling
./run_denovo_length.sh

# Get filters and annotations
./filter_and_annotate.sh
./apply_filters.sh

##### TODO change below to use filtered files
# Summarize families
./run_summarize_families2.sh

# Get all mutations
./get_mutations.sh
