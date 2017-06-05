#!/bin/bash

AWS_ACCESS_KEY=$(cat ~/.aws/credentials | grep "aws_access" | cut -f 3 -d' ')
AWS_SECRET_KEY=$(cat ~/.aws/credentials | grep "aws_secret" | cut -f 3 -d' ')
NUMPROC=3

for batch in $(ls -l batches/ | awk '{print $NF}' | grep -v bamfilesaa)
do
    bampaths=s3://ssc-psmc/batches/${batch}
    ./launch_aws.sh ${bampaths} ${AWS_ACCESS_KEY} ${AWS_SECRET_KEY} ${NUMPROC} micro_key
done

# Need to rerun:
# bamfilesbz
# bamfilesca
# bamfilescb
