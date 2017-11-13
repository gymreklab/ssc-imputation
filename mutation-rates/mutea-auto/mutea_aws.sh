#!/bin/bash

PARAMFILE=$1
BATCH=$2

source ${PARAMFILE}

die()
{
    BASE=$(basename -- "$0")
    echo "$BASE: error: $@" >&2
    exit 1
}

# Output paths
outdir=/mnt/batch_estimates 
outfile=ssc_hipstr_mutea_chrom${chrom}_batch${batchnum}.tab

# Check if outfile already on S3. If yes, don't recompute
x=$(aws s3 ls s3://ssc-mutea/batch_estimates/${outfile}.gz | awk '{print $NF}')
test -z $x || exit 0

# Download batch
aws s3 cp ${AWSBATCHPATH}/${BATCH} /mnt/batches/
batchpath=/mnt/batches/${BATCH}
batchnum=$(echo $BATCH | cut -f 2 -d'.')

# Download relevant VCF chunk
chrom=$(echo $BATCH | cut -f 1 -d'.')
start=$(head -n 1 $batchpath | cut -f 2)
end=$(tail -n 1 $batchpath | cut -f 3)
strvcf=/mnt/vcfs/${BATCH}.vcf.gz
tabix --print-header s3://ssc-strvcf/hipstr.chr${chrom}.asdt.vcf.gz ${chrom}:${start}-${end} | \
    bgzip -c > ${strvcf}
tabix -p vcf ${strvcf}

python /root/mutea-autosomal/mutea-auto/main_autosomal.py \
    --asdhet ${strvcf} --vcf \
    --out ${outdir}/${outfile} \
    --use-likelihoods --output-central-allele \
    --loci ${batchpath} \
    --min_samples ${MINSAMPLES} \
    --min_mu ${MINMU} --max_mu ${MAXMU} \
    --min_beta ${MINBETA} --max_beta ${MAXBETA} \
    --min_pgeom ${MINPGEOM} --max_pgeom ${MAXPGEOM} \
    --stderrs fisher || die "Mutea on ${BATCH} failed"

# Upload results to S3
gzip ${outdir}/${outfile}
aws s3 cp ${outdir}/${outfile}.gz s3://ssc-mutea/batch_estimates/${outfile}.gz

# Remove intermediate files so instance storage doesn't explode
rm ${strvcf}
rm ${batchpath}
rm ${outdir}/${outfile}.gz

exit 0