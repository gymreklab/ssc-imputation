#!/bin/bash

inline=$(cat -)

chrom=$(echo $inline | cut -f 1 -d' ')
pos=$(echo $inline | cut -f 2 -d' ')
period=$(echo $inline | cut -f 3 -d' ')
pnomut=$(echo $inline | cut -f 4 -d' ')
pmutaff=$(echo $inline | cut -f 5 -d' ')
pmutunaff=$(echo $inline | cut -f 6 -d' ')
mendlengths=$(echo $inline | cut -f 7 -d' ')
mendlengthsother=$(echo $inline | cut -f 8 -d' ')

affecteddata=$(echo $inline | cut -f 9 -d' ')
unaffecteddata=$(echo $inline | cut -f 10 -d' ')
motherdata=$(echo $inline | cut -f 11 -d' ')
fatherdata=$(echo $inline | cut -f 12 -d' ')

filter=$(echo $inline | cut -f 13 -d' ')
denovofilter=$(echo $inline | cut -f 14 -d' ')
poo=$(echo $inline | cut -f 15 -d' ')
newallele=$(echo $inline | cut -f 16 -d' ')
ref=$(echo $inline | cut -f 17 -d' ')
alt=$(echo $inline | cut -f 18 -d' ')

# Get alleles
allalleles=${ref},${alt}
m1=$(echo ${motherdata} | cut -d':' -f 1 | cut -d'|' -f 1 | awk '{print $0+1}')
m2=$(echo ${motherdata} | cut -d':' -f 1 | cut -d'|' -f 2 | awk '{print $0+1}')
m1allele=$(echo ${allalleles} | cut -d',' -f ${m1})
m2allele=$(echo ${allalleles} | cut -d',' -f ${m2})
f1=$(echo ${fatherdata} | cut -d':' -f 1 | cut -d'|' -f 1 | awk '{print $0+1}')
f2=$(echo ${fatherdata} | cut -d':' -f 1 | cut -d'|' -f 2 | awk '{print $0+1}')
f1allele=$(echo ${allalleles} | cut -d',' -f ${f1})
f2allele=$(echo ${allalleles} | cut -d',' -f ${f2})
a1=$(echo ${affecteddata} | cut -d':' -f 1 | cut -d'|' -f 1 | awk '{print $0+1}')
a2=$(echo ${affecteddata} | cut -d':' -f 1 | cut -d'|' -f 2 | awk '{print $0+1}')
a1allele=$(echo ${allalleles} | cut -d',' -f ${a1})
a2allele=$(echo ${allalleles} | cut -d',' -f ${a2})
u1=$(echo ${unaffecteddata} | cut -d':' -f 1 | cut -d'|' -f 1 | awk '{print $0+1}')
u2=$(echo ${unaffecteddata} | cut -d':' -f 1 | cut -d'|' -f 2 | awk '{print $0+1}')
u1allele=$(echo ${allalleles} | cut -d',' -f ${u1})
u2allele=$(echo ${allalleles} | cut -d',' -f ${u2})

echo "##### ${chrom}:${pos} - period=${period} - parentorigin=${poo}####"
echo "p_nomut=${pnomut}, p_mutaff=${pmutaff}, p_mutunaff=${pmutunaff}"
echo "mend lengths = ${mendlengths}, ${mendlengthsother}"
echo ""

echo "Mother"
echo "Allele1 - ${m1allele}"
echo "Allele2 - ${m2allele}"
echo "Genotype="$(echo ${motherdata} | cut -d':' -f 1) \
    "GB="$(echo ${motherdata} | cut -d':' -f 7) \
    "Q="$(echo ${motherdata} | cut -d':' -f 2) \
    "DSTUTTER="$(echo ${motherdata} | cut -d':' -f 3) \
    "DP="$(echo ${motherdata} | cut -d':' -f 4) \
    "DFLANKINDEL="$(echo ${motherdata} | cut -d':' -f 5)
echo "mallreads - mother"
mall=$(echo ${motherdata} | cut -d':' -f 6)
for i in $(seq 1 $(echo $mall | awk -F";" '{print NF}'))
do
    allele=$(echo $mall | cut -d';' -f ${i})
    echo "Allele "$(echo $allele | sed 's/|/: /')" reads"
done
echo "allreads - mother"
all=$(echo ${motherdata} | cut -d':' -f 8)
for i in $(seq 1 $(echo $all | awk -F";" '{print NF}'))
do
    allele=$(echo $all | cut -d';' -f ${i})
    echo "Allele "$(echo $allele | sed 's/|/: /')" reads"
done

echo
echo "Father"
echo "Allele1 - ${f1allele}"
echo "Allele2 - ${f2allele}"
echo "Genotype="$(echo ${fatherdata} | cut -d':' -f 1) \
    "GB="$(echo ${fatherdata} | cut -d':' -f 7) \
    "Q="$(echo ${fatherdata} | cut -d':' -f 2) \
    "DSTUTTER="$(echo ${fatherdata} | cut -d':' -f 3) \
    "DP="$(echo ${fatherdata} | cut -d':' -f 4) \
    "DFLANKINDEL="$(echo ${fatherdata} | cut -d':' -f 5)
echo "mallreads - father"
mall=$(echo ${fatherdata} | cut -d':' -f 6)
for i in $(seq 1 $(echo $mall | awk -F";" '{print NF}'))
do
    allele=$(echo $mall | cut -d';' -f ${i})
    echo "Allele "$(echo $allele | sed 's/|/: /')" reads"
done
echo "allreads - father"
all=$(echo ${fatherdata} | cut -d':' -f 8)
for i in $(seq 1 $(echo $all | awk -F";" '{print NF}'))
do
    allele=$(echo $all | cut -d';' -f ${i})
    echo "Allele "$(echo $allele | sed 's/|/: /')" reads"
done

echo
echo "Affected"
echo "Allele1 - ${a1allele}"
echo "Allele2 - ${a2allele}"
echo "Genotype="$(echo ${affecteddata} | cut -d':' -f 1) \
    "GB="$(echo ${affecteddata} | cut -d':' -f 7) \
    "Q="$(echo ${affecteddata} | cut -d':' -f 2) \
    "DSTUTTER="$(echo ${affecteddata} | cut -d':' -f 3) \
    "DP="$(echo ${affecteddata} | cut -d':' -f 4) \
    "DFLANKINDEL="$(echo ${affecteddata} | cut -d':' -f 5)
echo "mallreads - affected"
mall=$(echo ${affecteddata} | cut -d':' -f 6)
for i in $(seq 1 $(echo $mall | awk -F";" '{print NF}'))
do
    allele=$(echo $mall | cut -d';' -f ${i})
    echo "Allele "$(echo $allele | sed 's/|/: /')" reads"
done
echo "allreads - affected"
all=$(echo ${affecteddata} | cut -d':' -f 8)
for i in $(seq 1 $(echo $all | awk -F";" '{print NF}'))
do
    allele=$(echo $all | cut -d';' -f ${i})
    echo "Allele "$(echo $allele | sed 's/|/: /')" reads"
done

echo
echo "Unaffected"
echo "Allele1 - ${u1allele}"
echo "Allele2 - ${u2allele}"
echo "Genotype="$(echo ${unaffecteddata} | cut -d':' -f 1) \
    "GB="$(echo ${unaffecteddata} | cut -d':' -f 7) \
    "Q="$(echo ${unaffecteddata} | cut -d':' -f 2) \
    "DSTUTTER="$(echo ${unaffecteddata} | cut -d':' -f 3) \
    "DP="$(echo ${unaffecteddata} | cut -d':' -f 4) \
    "DFLANKINDEL="$(echo ${unaffecteddata} | cut -d':' -f 5)
echo "mallreads - unaffected"
mall=$(echo ${unaffecteddata} | cut -d':' -f 6)
for i in $(seq 1 $(echo $mall | awk -F";" '{print NF}'))
do
    allele=$(echo $mall | cut -d';' -f ${i})
    echo "Allele "$(echo $allele | sed 's/|/: /')" reads"
done
echo "allreads - unaffected"
all=$(echo ${unaffecteddata} | cut -d':' -f 8)
for i in $(seq 1 $(echo $all | awk -F";" '{print NF}'))
do
    allele=$(echo $all | cut -d';' -f ${i})
    echo "Allele "$(echo $allele | sed 's/|/: /')" reads"
done
