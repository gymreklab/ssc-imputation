# Get genome-wide file with features + observed rates
./get_features.sh # TODO
./combine_features_estimates.sh # TODO

# Get set of loci to train on
./get_training_loci.sh # TODO

# Build model and get per-locus constraint scores
runipy BuildConstraintModel.ipynb # TODO
