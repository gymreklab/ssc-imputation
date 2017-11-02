# Get genome-wide file with features + observed rates
./get_features.sh
./combine_features_estimates.sh

# Get set of loci to train on
./get_training_loci.sh

# Build model and get per-locus constraint scores
runipy BuildConstraintModel.ipynb # TODO
