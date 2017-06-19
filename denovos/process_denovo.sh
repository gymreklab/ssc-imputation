#!/bin/bash

source params.sh

chrom=$1
family=$2
father=$3
mother=$4
affected=$5
unaffected=$6

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
}

die()
{
    BASE=$(basename -- "$0")
    echo "$BASE: error: $@" >&2
    exit 1
}

echo "Family" ${family}
echo "Father" ${father}
echo "Mother" ${mother}
echo "Affected" ${affected}
echo "Unaffected" ${unaffected}

test -z ${family} && usage
test -z ${father} && usage
test -z ${mother} && usage
test -z ${affected} && usage
test -z ${unaffected} && usage
