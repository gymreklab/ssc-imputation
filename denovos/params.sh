DENOVOFINDER=/home/mgymrek/workspace/HipSTR-mgymrek/DenovoFinder
FAMFILE=/home/mgymrek/workspace/ssc-imputation/denovos/pedigree.fam
PARENTS=/home/mgymrek/workspace/ssc-imputation/metadata/ssc_parent_ids.txt
WORKDIR=/home/mgymrek/workspace/ssc-imputation/denovos
TMPDIR=/storage/mgymrek/ssc-denovos/test/tmp
BASEOUTDIR=/storage/mgymrek/ssc-denovos/denovos/
startchrom=21 # TODO change
endchrom=21 # TODO change

# Denovofinder options
WINDOWSIZE=50000
MAXALLELES=25

# Filtering options for individual genotypes 
MINQ=0.9 # Minimum genotype quality score for all family members
MINST=1.0 # Filter if %sttuter/indel reads > this
MINSPAN=10 # Require at least this many spanning reads
MINSUPP=0.15 # Each allele must be supported by this percent of reads
PARCOUNT=1000 # If parent has > this many reads supporing new allele, filter

# Summary filter levels
PROBTHRESH=0.9 # require max posterior prob to be at least this
MINFAM=100 # Require this many families to have data
