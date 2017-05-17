#!/bin/bash

set -e

BAMPATHS=$1
AWS_ACCESS_KEY=$2
AWS_SECRET_KEY=$3
NUMPROC=$4
KEYNAME=$5

usage()
{
    BASE=$(basename -- "$0")
    echo "Launch amazon instance to get diploid consensus sequence for each ssc bam parent
Usage:
    $BASE <bampaths> <aws_access_key> <aws_secret_key> <numproc> <keyname>
"
    exit 1
}
test -z ${BAMPATHS} && usage
test -z ${AWS_ACCESS_KEY} && usage
test -z ${AWS_SECRET_KEY} && usage
test -z ${NUMPROC} && usage
test -z ${KEYNAME} && usage

# Instance details
SPOT_PRICE=0.50
INSTANCE_TYPE=c4.large
IMAGE_ID=ami-80861296

STARTUP_SCRIPT=$(cat /home/mgymrek/workspace/ssc-imputation/psmc/run_from_aws.sh | \
    sed "s/\$1/${AWS_ACCESS_KEY}/" | sed "s~\$2~${AWS_SECRET_KEY}~" | \
    sed "s~\$3~${BAMPATHS}~" | \
    sed "s/\$4/${NUMPROC}/")
STARTUP_SCRIPT_ENCODE="$(echo "${STARTUP_SCRIPT}" | base64 -w 0)"

LAUNCH_SPEC="{\"ImageId\":\"${IMAGE_ID}\",\"SecurityGroupIds\":[\"sg-5e914222\"], \"KeyName\":\"${KEYNAME}\",\"InstanceType\":\"${INSTANCE_TYPE}\", \"UserData\":\"${STARTUP_SCRIPT_ENCODE}\"}"

aws ec2 request-spot-instances \
    --spot-price ${SPOT_PRICE} \
    --instance-count 1 \
    --availability-zone-group us-east-1a \
    --type one-time \
    --launch-specification "${LAUNCH_SPEC}"

