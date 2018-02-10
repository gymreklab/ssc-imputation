OUTDIR=/scratch/gymrek/genomewide-imputation-pgc/
SCRATCHDIR=${OUTDIR}/tmp/
COVARFILE=/home/pgcdac/DWFV2CJb8Piv_0116_pgc_data/scz/wave2/v1/prune.bfile.cobg.PGC_SCZ49.sh2.menv.mds_cov
COVARCOLS=C1,C2,C3,C4,C5,C6,C7,C9,C15,C18 # based on recommendation from S. Ripke

COHORTFILES=/home/gymrek/pgc_imputation/pgc_eur_mergelist.txt

BEAGLE=/home/gymrek/bin/beagle.08Jun17.d8b.jar
HIPREF=/home/gymrek/dbase/GRCh37.hipstr_reference.bed
REGIONS=pgc_gwas_regions.bed
WINDOW=50000