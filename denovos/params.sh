DENOVOFINDER=/home/mgymrek/workspace/HipSTR-mgymrek/DenovoFinder
FAMFILE=/home/mgymrek/workspace/ssc-imputation/denovos/pedigree.fam
WORKDIR=/home/mgymrek/workspace/ssc-imputation/denovos
TMPDIR=/storage/mgymrek/ssc-denovos/test/tmp
OUTDIR=/storage/mgymrek/ssc-denovos/test/ # TODO change
startchrom=12 # TODO change
endchrom=12 # TODO change

# Filtering options
MINQ=0.9 # Minimum genotype quality score for all family members
MINST=0.2 # Filter if %sttuter/indel reads > this
MINSPAN=10 # Require at least this many spanning reads
MINSUPP=0.2 # Each allele must be supported by this percent of reads
PARCOUNT=3 # If parent has > this many reads supporing new allele, filter
PROBTHRESH=0.9 # require max posterior prob to be at least this
