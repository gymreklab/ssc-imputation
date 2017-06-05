#!/bin/bash

LOGDIR=$1

# Get broken ones
grep -A 1 "bgzf_read_block error" ${LOGDIR}/*.log | grep upload | cut -f 2 -d':' | sed 's/ .\///' | sed 's/ to s3//' > broken.txt
