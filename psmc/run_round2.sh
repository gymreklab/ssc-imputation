#!/bin/bash

LAB_AWS_ACCESS_KEY=$(cat ~/.aws/credentials | grep "aws_access" | head -n 1 | cut -f 3 -d' ')
LAB_AWS_SECRET_KEY=$(cat ~/.aws/credentials | grep "aws_secret" | head -n 1 | cut -f 3 -d' ')
SSC_AWS_ACCESS_KEY=$(cat ~/.aws/ssc_credentials | grep "aws_access" | cut -f 3 -d' ')
SSC_AWS_SECRET_KEY=$(cat ~/.aws/ssc_credentials | grep "aws_secret" | cut -f 3 -d' ')

NUMPROC=3

for batch in $(ls -l batches_round2/ | awk '{print $NF}' | grep bamfiles2 | grep -v bamfiles2aa)
do
    bampaths=s3://ssc-psmc/batches_round2/${batch}
    ./launch_aws_round2.sh ${bampaths} \
	${LAB_AWS_ACCESS_KEY} ${LAB_AWS_SECRET_KEY} \
	${SSC_AWS_ACCESS_KEY} ${SSC_AWS_SECRET_KEY} \
	${NUMPROC} micro_key
done
