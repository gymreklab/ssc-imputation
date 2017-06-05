#!/bin/bash

LOGDIR=$1

# Get broken ones
cat ${LOGDIR}/*.log | grep -A 1 -B 1 "<mpileup>" | grep upload | cut -f 2 -d':' | sed 's/ .\///' | sed 's/ to s3//' > done.txt
