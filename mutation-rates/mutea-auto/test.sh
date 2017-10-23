#!/bin/bash

#SBATCH -A ddp268
#SBATCH -t 100
#SBATCH --mem=5G
#SBATCH -p shared
#SBATCH --job-name=testmutea
#SBATCH --get-user-env
#SBATCH -o testmutea.out
#SBATCH -e testmutea.err

time ./run_mutea_autosomal.sh 2 18 "--maxloci 10"
