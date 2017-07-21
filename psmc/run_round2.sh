#!/bin/bash

BATCHDIR=batches_round3
BATCHSIZE=160
DONE=all_completed_jobs.txt

# Get completed jobs
#cat completed_round1.txt > ${DONE}
#aws s3 ls s3://ssc-psmc/round2/ | grep "final.bam" | awk '{print $NF}' >> ${DONE}

# Get remaining jobs, split into batches
#./get_batches_round2.sh ${BATCHDIR} ${BATCHSIZE} ${DONE}

LAB_AWS_ACCESS_KEY=$(cat ~/.aws/credentials | grep "aws_access" | head -n 1 | cut -f 3 -d' ')
LAB_AWS_SECRET_KEY=$(cat ~/.aws/credentials | grep "aws_secret" | head -n 1 | cut -f 3 -d' ')
SSC_AWS_ACCESS_KEY=$(cat ~/.aws/ssc_credentials | grep "aws_access" | cut -f 3 -d' ')
SSC_AWS_SECRET_KEY=$(cat ~/.aws/ssc_credentials | grep "aws_secret" | cut -f 3 -d' ')

NUMPROC=3

for batch in $(ls -l ${BATCHDIR}/ | awk '{print $NF}' | grep bamfiles2)
do
    bampaths=s3://ssc-psmc/${BATCHDIR}/${batch}
    ./launch_aws_round2.sh ${bampaths} \
	${LAB_AWS_ACCESS_KEY} ${LAB_AWS_SECRET_KEY} \
	${SSC_AWS_ACCESS_KEY} ${SSC_AWS_SECRET_KEY} \
	${NUMPROC} micro_key
done
