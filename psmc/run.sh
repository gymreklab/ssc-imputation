#!/bin/bash

BAMPATHS=s3://ssc-psmc/batches/bamfilesaa
AWS_ACCESS_KEY=$(cat ~/.aws/credentials | grep "aws_access" | cut -f 3 -d' ')
AWS_SECRET_KEY=$(cat ~/.aws/credentials | grep "aws_secret" | cut -f 3 -d' ')
NUMPROC=3

./launch_aws.sh ${BAMPATHS} ${AWS_ACCESS_KEY} ${AWS_SECRET_KEY} ${NUMPROC} micro_key
