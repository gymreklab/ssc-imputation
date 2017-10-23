#!/bin/bash

#SBATCH -A ddp268
#SBATCH -t 100
#SBATCH --mem=5G
#SBATCH -p shared
#SBATCH --job-name=testmutea
#SBATCH --get-user-env
#SBATCH -o testcodis.out
#SBATCH -e testcodis.err

source params.sh

python ${MUTEADIR}/mutea-auto/main_autosomal.py \
    --asdhet ${ASDT_VCF} --vcf --asdhetdir \
    --out ${OUTDIR}/test/ssc_hipstr_mutea_codis.tab \
    --usestutter ${OUTDIR}/test/stutter_codis.bed \
    --loci ${CODIS} \
    --min_samples ${MINSAMPLES} \
    --min_mu ${MINMU} --max_mu ${MAXMU} \
    --min_beta ${MINBETA} --max_beta ${MAXBETA} \
    --min_pgeom ${MINPGEOM} --max_pgeom ${MAXPGEOM} \
    --stderrs fisher 
