#!/bin/bash

set -e

BAMPATHS=$1
LAB_AWS_ACCESS_KEY=$2
LAB_AWS_SECRET_KEY=$3
SSC_AWS_ACCESS_KEY=$4
SSC_AWS_SECRET_KEY=$5
NUMPROC=$6
KEYNAME=$7

usage()
{
    BASE=$(basename -- "$0")
    echo "Launch amazon instance to get diploid consensus sequence for each ssc bam parent
Usage:
    $BASE <bampaths> <lab_aws_access_key> <lab_aws_secret_key> <ssc_aws_access_key> <ssc_aws_secret_key> <numproc> <keyname>
"
    exit 1
}
test -z ${BAMPATHS} && usage
test -z ${LAB_AWS_ACCESS_KEY} && usage
test -z ${LAB_AWS_SECRET_KEY} && usage
test -z ${SSC_AWS_ACCESS_KEY} && usage
test -z ${SSC_AWS_SECRET_KEY} && usage
test -z ${NUMPROC} && usage
test -z ${KEYNAME} && usage

# Instance details
SPOT_PRICE=0.05
INSTANCE_TYPE=c4.large
IMAGE_ID=ami-80861296

STARTUP_SCRIPT=$(cat /home/mgymrek/workspace/ssc-imputation/psmc/run_from_aws_round2.sh | \
    sed "s/\$1/${LAB_AWS_ACCESS_KEY}/" | sed "s~\$2~${LAB_AWS_SECRET_KEY}~" | \
    sed "s/\$3/${SSC_AWS_ACCESS_KEY}/" | sed "s~\$4~${SSC_AWS_SECRET_KEY}~" | \
    sed "s~\$5~${BAMPATHS}~" | \
    sed "s/\$6/${NUMPROC}/")
STARTUP_SCRIPT_ENCODE="$(echo "${STARTUP_SCRIPT}" | base64 -w 0)"

LAUNCH_SPEC="{\"ImageId\":\"${IMAGE_ID}\",\"Placement\":{\"AvailabilityZone\": \"us-east-1b\"},\"SecurityGroupIds\":[\"sg-5e914222\"], \"KeyName\":\"${KEYNAME}\",\"InstanceType\":\"${INSTANCE_TYPE}\", \"UserData\":\"${STARTUP_SCRIPT_ENCODE}\"}"

aws ec2 request-spot-instances \
    --spot-price ${SPOT_PRICE} \
    --instance-count 1 \
    --type one-time \
    --launch-specification "${LAUNCH_SPEC}"

