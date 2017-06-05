#!/bin/bash

LOGDIR=$1

aws ec2 describe-instances --filter "Name=instance-type,Values=c4.large" --out json | grep PublicDnsName | cut -f 2 -d ':' | sed 's/ \"//' | sed 's/\",//' | uniq > ${LOGDIR}/ec2ids.txt

for eid in $(cat ${LOGDIR}/ec2ids.txt)
do
    scp -o StrictHostKeyChecking=no -i ~/keys/micro_key.pem ubuntu@${eid}:/var/log/cloud-init-output.log ${LOGDIR}/${eid}.log 
done

