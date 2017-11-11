#!/bin/bash

source params.sh

BATCHDIR=${OUTDIR}/batches/

aws s3 sync ${BATCHDIR} s3://ssc-mutea/batches/
