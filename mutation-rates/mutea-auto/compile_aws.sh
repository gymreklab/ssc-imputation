#!/bin/bash

#SBATCH -A ddp268
#SBATCH -p shared
#SBATCH -t 100
#SBATCH --get-user-env

source params.sh

# Download all batch estimates from AWS
aws s3 sync s3://ssc-mutea/batch_estimates/ ${OUTDIR}/batch_estimate_aws/

# Gather completed
zcat ${OUTDIR}/batch_estimate_aws/*.gz | \
    sort -k1,1 -k2,2n | uniq | bgzip -c > ${OUTDIR}/genome-wide/ssc_mutea_auto_unfiltered_aws.bed.gz
tabix -f -p bed ${OUTDIR}/genome-wide/ssc_mutea_auto_unfiltered_aws.bed.gz
