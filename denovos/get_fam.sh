#!/bin/bash

IDFILE=~/workspace/ssc-imputation/metadata/ssc_family_ids.txt
PILOTFILE=pilot_families.txt

# Get fam from SSC ids file
families=$(cat ${IDFILE} | cut -f 1 -d'.' | grep -v -w -f ${PILOTFILE} | sort | uniq)
for family in ${families}
do
    father=$(cat ${IDFILE} | grep ${family}.fa | cut -f 2)
    mother=$(cat ${IDFILE} | grep ${family}.mo | cut -f 2)
    affected=$(cat ${IDFILE} | grep ${family}.p1 | cut -f 2)
    unaffected=$(cat ${IDFILE} | grep ${family}.s1 | cut -f 2)
    # Write .fam file
    echo $family $father 0 0 1 -9
    echo $family $mother 0 0 2 -9
    echo $family $affected $father $mother 0 2
    echo $family $unaffected $father $mother 0 1
done
