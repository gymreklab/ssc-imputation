#!/bin/bash

# TODO first download all s3 logs

# Get broken ones
grep -A 1 "bgzf_read_block error" *.log | grep upload | cut -f 2 -d':' | sed 's/ .\///' | sed 's/ to s3//' > broken.txt
