#!/bin/bash

# Example: ./process_denovo.sh 12 14652 SSC12871 SSC12867 SSC12865 SSC12872


source params.sh
cd ${WORKDIR}

usage()
{
    BASE=$(basename -- "$0")
    echo "Process de novo calls from nuclear family
Usage:
    $BASE <chrom> <family> <father> <mother> <affected> <unaffected>
Get all calls for the family in tab format. Downstream scripts will filter
"
exit 1
}

die()
{
    BASE=$(basename -- "$0")
    echo "$BASE: error: $@" >&2
    exit 1
}

chrom=$1
family=$2
father=$3
mother=$4
affected=$5
unaffected=$6

test -z ${chrom} && usage
test -z ${family} && usage
test -z ${father} && usage
test -z ${mother} && usage
test -z ${affected} && usage
test -z ${unaffected} && usage

echo "Chrom" ${chrom}
echo "Family" ${family}
echo "Father" ${father}
echo "Mother" ${mother}
echo "Affected" ${affected}
echo "Unaffected" ${unaffected}

outfile=${BASEOUTDIR}/bychrom/SSC_denovo_${family}_${chrom}.tab
denovovcf=${BASEOUTDIR}/denovofinder_trio/denovofinder_chr${chrom}.vcf.gz
gtvcf=/storage/s1saini/hipstr_genomewide/chr${chrom}/hipstr_calls_${chrom}.vcf.gz

echo ${TMPDIR}
# Get denovo info - NOMUT:DENOVO:OTHER
echo "[process_denovo_trio.sh] Getting de novo info... "
bcftools query -s${unaffected},${affected} -f '%POS[\t%NOMUT\t%DENOVO\t%OTHER]\n' \
    ${denovovcf} \
    > ${TMPDIR}/${family}.${chrom}.denovo || \
    die "Could not get denovo info ${family}"

# Get genotype info
echo "[process_denovo_trio.sh] Getting genotype info... "
bcftools query -s${unaffected},${affected},${father},${mother} \
    -f '%POS\t%PERIOD\t%REF\t%ALT[\t%GT:%Q:%DSTUTTER:%DP:%DFLANKINDEL:%MALLREADS:%GB:%ALLREADS]\n' \
    ${gtvcf} \
    > ${TMPDIR}/${family}.${chrom}.gts || \
    die "Could not get gt info ${family}"

# Get output file
echo "[process_denovo_trio.sh] Getting output file..."
echo "chrom,pos,nomut_ll_unaffected,denovo_ll_unaffected,other_ll_unaffected,nomut_ll_affected,denovo_ll_affected,other_ll_affected,period,ref,alt,gtdata_unaffected,gtdata_affected,gtdata_father,gtdata_mother" | \
    sed 's/,/\t/g' > ${outfile}
join -t $'\t' ${TMPDIR}/${family}.${chrom}.denovo ${TMPDIR}/${family}.${chrom}.gts | \
    awk -v"chrom=$chrom" '{print chrom "\t" $0}' >> ${outfile}

exit 0
