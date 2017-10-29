# Split loci and stutter to batches
./get_batches.sh

# Run mutea autosomal on all batches
./run_mutea_all.sh

# Calibrate stderrs using CODIS
./test_codis.sh
./calibrate_errors.sh

# Gather batches
./gather_completed.sh

# Filter and scale
./filter_and_scale.sh

# Summarize mutation models
./summarize_mutation_models.sh

###################
