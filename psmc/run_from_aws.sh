#!/bin/bash

AWS_ACCESS_KEY=$1
AWS_SECRET_KEY=$2
BAMPATHS=$3
NUMPROC=$4

AWS_DIR=${HOME}/.aws
AWS_CONFIG_FILE=${AWS_DIR}/config
AWS_CRED_FILE=${AWS_DIR}/credentials
JOBSFILE=jobs.txt

OUTBUCKET=s3://ssc-psmc
REFFA=/mnt/tmp/Homo_sapiens_assembly19.fasta

usage()
{
    BASE=$(basename -- "$0")
    echo "Run PSMC analysis
Usage:
    $BASE <aws access key> <aws secret key> <bamfileslist> <numproc>
       - aws access key and aws secret keys are for AWS configuration
       - bamfileslist has list of S3 paths to bams to process. It must be stored on S3
       - numproc gives number of jobs to run at once
Does the following:
1. Set up AWS configuration
2. Download necessary files
3. Create jobs
4. Run those jobs
5. Upload results to S3 bucket
6. Terminate
"
    terminate
    exit 1
}

terminate() {
    INSTANCE_ID=$(curl http://169.254.169.254/latest/meta-data/instance-id)
    # Terminate instance
    echo "Terminating instance ${INSTANCE_ID}"
    aws ec2 terminate-instances --output table --instance-ids ${INSTANCE_ID}
    exit 1 # shouldn't happen
}

test -z ${AWS_ACCESS_KEY} && usage
test -z ${AWS_SECRET_KEY} && usage
test -z ${BAMPATHS} && usage
test -z ${NUMPROC} && usage

die()
{
    BASE=$(basename -- "$0")
    echo "$BASE: error: $@" >&2
    terminate
    exit 1
}

# Install things
apt-get update || die "Could not update"
apt-get -y install awscli || die "Could not install aws"
apt-get -y install git || die "Could not install git"

# Set up AWS credentials
mkdir -p ${AWS_DIR} || die "Could not create AWS dir"
echo "[default]" > ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CONFIG_FILE}"
echo "output = table" >> ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CONFIG_FILE}"
echo "region = us-east-1" >> ${AWS_CONFIG_FILE}  || die "Could not write to ${AWS_CONFIG_FILE}"
echo "aws_secret_access_key = ${AWS_SECRET_KEY}" >> ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "aws_access_key_id = ${AWS_ACCESS_KEY}" >> ${AWS_CONFIG_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "[default]" > ${AWS_CRED_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "aws_secret_access_key = ${AWS_SECRET_KEY}" >> ${AWS_CRED_FILE} || die "Could not write to ${AWS_CRED_FILE}"
echo "aws_access_key_id = ${AWS_ACCESS_KEY}" >> ${AWS_CRED_FILE} || die "Could not write to ${AWS_CRED_FILE}"

# Get github
cd ~
git clone https://github.com/gymreklab/ssc-imputation || die "Could not clone github repo"

# Download files
mkdir /mnt/tmp || die "Could not make tmp directory"
aws s3 cp ${OUTBUCKET}/Homo_sapiens_assembly19.fasta ${REFFA}
mkdir /mnt/tmp/consensus/
aws s3 cp ${BAMPATHS} /mnt/tmp/bamfiles.txt

# Create jobs
cd /mnt/tmp/consensus/
for bamfile in $(cat /mnt/tmp/bamfiles.txt)
do
    outfile=/mnt/tmp/consensus/$(basename ${bam}).fq.gz
    echo "~/workspace/ssc-imputation/psmc/get_consensus.sh ${bamfile} ${OUTBUCKET} ${outfile} ${REFFA}" >> ${JOBSFILE}
done

# Run jobs
cat ${JOBSFILE} | xargs -n 1 -P${NUMPROC} -I% bash -c %

terminate

