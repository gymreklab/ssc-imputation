#!/bin/bash

source params.sh

cat ${DIR}/snp_str_ld.tab | datamash -g 2 max 5 -f