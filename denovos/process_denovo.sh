#!/bin/bash

# Example: ./process_denovo.sh 12 14652 SSC12871 SSC12867 SSC12865 SSC12872


source params.sh
cd ${WORKDIR}

usage()
{
    BASE=$(basename -- "$0")
    echo "Process de novo calls from single SSC family
Usage:
    $BASE <chrom> <family> <father> <mother> <affected> <unaffected>
Does the following:
1. Get all calls for each family in tab format
2. Apply python script to perform basic filters on calls
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

outfile=${OUTDIR}/bychrom/SSC_denovo_${family}_${chrom}.tab
denovovcf=${OUTDIR}/test_denovos_withsnps.vcf.gz # TODO change
gtvcf=/storage/mgymrek/ssc-denovos/test/hipstr_calls_12_145.vcf.gz # TODO change

# Get denovo info CHILDREN:NOMUT:ANYMUT:DENOVO:OTHER
echo "[process_denovo.sh] Getting de novo info..."
bcftools query -s${family} -f '[%POS\t%NOMUT\t%ANYMUT\t%DENOVO\t%OTHER\n]' \
    ${denovovcf} | sort -k1,1 \
    > ${TMPDIR}/${family}.${chrom}.denovo || die "Could not get denovo info ${family}"

# Get genotype info
echo "[process_denovo.sh] Getting genotype info..."
bcftools query -s${affected},${unaffected},${father},${mother} \
    -f '%POS\t%PERIOD\t%REF\t%ALT[\t%GT:%Q:%DSTUTTER:%DP:%DFLANKINDEL:%MALLREADS:%GB:%ALLREADS]\n' \
    ${gtvcf} | sort -k1,1 \
    > ${TMPDIR}/${family}.${chrom}.gts || die "Could not get gt info ${family}"

# Get output file
echo "[process_denovo.sh] Getting output file..."
echo "chrom,pos,nomut_ll,anymut_ll,denovo_ll,other_ll,period,ref,alt,gtdata_affected,gtdata_unaffected,gtdata_father,gtdata_mother" | \
    sed 's/,/\t/g' > ${outfile}
join -t $'\t' ${TMPDIR}/${family}.${chrom}.denovo ${TMPDIR}/${family}.${chrom}.gts | \
    awk -v"chrom=$chrom" '{print chrom "\t" $0}' >> ${outfile}

echo "[process_denovo.sh] Running filter_denovos.py..."
./filter_denovos.py \
    --infile ${outfile} \
    --outfile ${outfile}.filtered \
    --minq-child ${MINQ} --minq-father ${MINQ} --minq-mother ${MINQ} \
    --minst-child ${MINST} --minst-father ${MINST} --minst-mother ${MINST} \
    --min-spanning-child ${MINSPAN} --min-spanning-father ${MINSPAN} --min-spanning-mother ${MINSPAN} \
    --min-percsupp-child ${MINSUPP} --min-percsupp-father ${MINSUPP} --min-percsupp-mother ${MINSUPP} \
    --parent-allele-count ${PARCOUNT} \
    --unit

exit 0
