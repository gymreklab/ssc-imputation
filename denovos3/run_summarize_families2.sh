#!/bin/bash
#SBATCH -A csd568
#SBATCH --get-user-env 
#SBATCH --job-name Summarize_SSC_Families
#SBATCH -p shared
#SBATCH -t 1000
#SBATCH --mem=16gb
#SBATCH -o /oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered_mutations/logs/Summarize_SSC_Families.log
#SBATCH -e /oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered_mutations/logs/Summarize_SSC_Families.log                                                          
#SBATCH --mail-user=ilmitra@ucsd.edu
#SBATCH --mail-type=ALL

date;hostname;pwd

source /home/ileena/ssc-imputation/denovos3/params.sh

AFTERFILTER=1

#AFTERFILTER=$1
if [ x$AFTERFILTER == x0 ]
then
    DIR=$OUTDIR
else
    DIR=$FINALOUTDIR
fi

mkdir $DIR/tmp

for period in $(seq 1 6)
do
	echo "/home/ileena/ssc-imputation/denovos3/summarize_families2.sh ${period}" $AFTERFILTER
done | xargs -P6 -I% -n1 sh -c "%"

cat ${DIR}/denovos_bylength_bychild2_period* | \
    head -n 1 > ${DIR}/denovos_bylength_bychild2.tab
cat ${DIR}/denovos_bylength_bychild2_period* | \
    awk '{print $NF "\t" $0}' | sort | \
    awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8}' | \
    grep -v family | datamash -g 1,2,3,4 sum 5 sum 6 sum 7 mean 8 | \
    cut -f 1 --complement >> ${DIR}/denovos_bylength_bychild2.tab
