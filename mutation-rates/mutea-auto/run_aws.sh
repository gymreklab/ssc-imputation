#!/bin/bash

source params.sh

AWS_ACCESS_KEY=$(cat ~/.aws/credentials | grep "aws_access" | cut -f 3 -d' ')
AWS_SECRET_KEY=$(cat ~/.aws/credentials | grep "aws_secret" | cut -f 3 -d' ')

for batch in $(ls -l ${OUTDIR}/superbatches_aws/ | '{print $NF}')
do
    batchpath=s3://ssc-mutea/superbatches/${batch}
    echo ./launch_aws.sh ${batchpath} ${AWS_ACCESS_KEY} ${AWS_SECRET_KEY} micro_key # TODO remove echo
done