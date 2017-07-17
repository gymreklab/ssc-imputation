#!/bin/bash

AWS_ACCESS_KEY=$(cat ~/.aws/ssc_credentials | grep "aws_access" | cut -f 3 -d' ')
AWS_SECRET_KEY=$(cat ~/.aws/ssc_credentials | grep "aws_secret" | cut -f 3 -d' ')
NUMPROC=3

for batch in $(ls -l batches_round2/ | awk '{print $NF}')
do
    bampaths=s3://ssc-psmc/batches_round2/${batch}
    ./launch_aws_round2.sh ${bampaths} ${AWS_ACCESS_KEY} ${AWS_SECRET_KEY} ${NUMPROC} micro_key
done
