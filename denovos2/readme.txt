# Run denovo calling
./run_denovo_length.sh

# Get filters and annotations
./filter_and_annotate.sh
./apply_filters.sh
./summarize_loci.sh

# Summarize families
./run_summarize_families2.sh 1

####### For intermediate analysis ####
./run_summarize_families2.sh 0 # Use to determine which posterior threshold to use

