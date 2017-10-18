# Get file with local TMRCA for each STR, scaled to genome-wide average of 21,000 generations
# Want format: chrom, start, end, tmrca, sample. One sorted and indexed file per chromosome
./get_str_tmrcas.sh

# Add TMRCA field to each VCF - TODO annotate_vcf.py
sbatch -A ddp268 -t 2000 --array=1-22 --job-name=annotatevcfs --get-user-env ./annotate_vcfs.sh
