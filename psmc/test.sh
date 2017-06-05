#!/bin/bash

OUTBUCKET=s3://ssc-psmc
OUTFILE=$1

x=$(aws s3 ls ${OUTBUCKET}/$(basename ${OUTFILE}) | awk '{print $NF}')
test -z $x || exit 0

echo "file not on aws. proceeding"
