#!/bin/bash

BAMFILES=../metadata/ssc_parent_bampaths.txt
split -l 20 ${BAMFILES} batches/bamfiles
aws s3 sync batches/ s3://ssc-psmc/batches/
