#!/bin/bash

AWS_ACCESS_KEY=$1
AWS_SECRET_KEY=$2
SUPERBATCHPATH=$3

HOMEDIR=/root/

AWS_DIR=${HOMEDIR}/.aws
AWS_CONFIG_FILE=${AWS_DIR}/config
AWS_CRED_FILE=${AWS_DIR}/credentials

OUTBUCKET=s3://ssc-mutea/batch_estimates

superbatch=$(basename $SUPERBATCHPATH)

usage()
{
    BASE=$(basename -- "$0")
    echo "Run MUTEA on SSC
Usage:
    $BASE <aws access key> <aws secret key> <superbatch>

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
    # Get log
    sudo aws s3 cp --output table /var/log/cloud-init-output.log ${OUTBUCKET}/log/${superbatch}.log
    # Terminate instance
    echo "Terminating instance ${INSTANCE_ID}"
#   sudo aws ec2 terminate-instances --output table --instance-ids ${INSTANCE_ID} #TODO uncomment
    exit 1 # shouldn't happen
}

test -z ${AWS_ACCESS_KEY} && usage
test -z ${AWS_SECRET_KEY} && usage
test -z ${SUPERBATCHPATH} && usage

die()
{
    BASE=$(basename -- "$0")
    echo "$BASE: error: $@" >&2
    terminate
    exit 1
}

# Install things
sudo apt-get update || die "Could not update"
sudo apt-get -y install awscli || die "Could not install aws"
sudo apt-get -y install git || die "Could not install git"
sudo apt-get -y install make gcc libz-dev libncurses5-dev libbz2-dev liblzma-dev libcurl3-dev libssl-dev autoconf || die "Could not install devtools"
sudo apt-get install -y python-setuptools python-dev build-essential || die "Could not install python"
sudo easy_install pip || die "Could not install pip"
sudo pip install joblib numpy scipy pytabix pyvcf statsmodels pycallgraph pysam || die "Could not install python libraries"

mkdir -p ${HOMEDIR}/source || die "Could not create source dir"
cd ${HOMEDIR}/source || die "Could not go to source dir"

cd ${HOMEDIR}/source
wget https://github.com/samtools/bcftools/releases/download/1.4.1/bcftools-1.4.1.tar.bz2
tar -xvf bcftools-1.4.1.tar.bz2
cd bcftools-1.4.1
make
sudo make install

cd ${HOMEDIR}
git clone https://github.com/samtools/htslib
cd htslib
autoheader
autoconf
./configure --enable-libcurl
make
sudo make install

cd ${HOMEDIR}
git clone https://github.com/samtools/samtools
cd samtools
autoconf -Wno-syntax 
./configure
make
sudo make install

# Set up AWS credentials
echo "Setting up AWS credentials in ${AWS_DIR}"
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
cd ${HOMEDIR}
git clone https://github.com/gymreklab/ssc-imputation || die "Could not clone github repo ssc-imputation"
cd ${HOMEDIR}
git clone https://github.com/gymreklab/mutea-autosomal || die "Could not clone github repo mutea-autosomal"

# Download files
sudo mkdir -p /mnt/tmp || die "Could not make tmp directory"
aws s3 cp ${SUPERBATCHPATH} /mnt/tmp/superbatch.txt

# Make directory for inputs/outputs
sudo mkdir -p /mnt/vcfs || die "Could not make vcfs directory"
sudo mkdir -p /mnt/batch_estimates || die "Could not make output directory"
sudo mkdir -p /mnt/batches/

# Run each job
source ${HOMEDIR}/ssc-imputation/mutation-rates/mutea-auto/params.sh
for batch in $(cat /mnt/tmp/superbatch.txt)
do
    # Download files from S3
    sudo ${HOMEDIR}/ssc-imputation/mutation-rates/mutea-auto/mutea_aws.sh \
	${HOMEDIR}/ssc-imputation/mutation-rates/mutea-auto/params.sh \
	${batch}
done

terminate

