# Get file with local TMRCA for each STR, scaled to genome-wide average of 21,000 generations
# Want format: chrom, start, end, tmrca, sample. One sorted and indexed file per chromosome
sbatch ./get_str_tmrcas_preprocess.sh
sbatch -p shared -A ddp268 -t 2800 --array=1-22 --job-name=getstrtmrca --get-user-env ./get_str_tmrcas.sh

# Add TMRCA field to each VCF
sbatch -p shared -A ddp268 -t 2800 --array=1-22 --job-name=annotatevcfs --get-user-env ./annotate_vcfs.sh
