#!/bin/bash

PICARD=/storage/resources/source/picard.jar
REF=/storage/resources/dbase/human/hs37d5/hs37d5

java -jar ${PICARD} CreateSequenceDictionary \
    R=${REF}.fa \
    O=${REF}.dict
