#!/bin/bash

# TODO first download all s3 logs

# Get broken ones
cat *.log | grep -A 1 -B 1 "<mpileup>" | grep upload | cut -f 2 -d':' | sed 's/ .\///' | sed 's/ to s3//' > done.txt
