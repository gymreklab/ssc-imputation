ASDT_VCF=/oasis/scratch/comet/mgymrek/temp_project/asdt_vcfs
OUTDIR=/oasis/projects/nsf/ddp268/mgymrek/ssc-quads/mutea-autosomal/
CONSTRAINTDIR=/oasis/projects/nsf/ddp268/mgymrek/ssc-quads/constraint/
LOGDIR=/oasis/projects/nsf/ddp268/mgymrek/ssc-quads/mutea-autosomal/log/
MUTEADIR=/home/mgymrek/workspace/mutea-autosomal
HIPREF=/oasis/projects/nsf/ddp268/mgymrek/dbase/GRCh37.hipstr_reference_sorted.bed
STUTTERFILE=/oasis/projects/nsf/ddp268/mgymrek/ssc-quads/hipstr-call-info/ssc_hipstr_stutter_hg19.bed
TMPLOC=/tmp

CODIS=/home/mgymrek/workspace/ssc-imputation/mutation-rates/mutea-auto/CODIS.bed
BATCHSIZE=5000
AWSBATCHSIZE=500 # size of individual batches
AWSSUPERBATCHSIZE=70 # each instance gets this many batches to process
AWSBATCHPATH=s3://ssc-mutea/batches

MINSAMPLES=50
MINMU=0.00000001
MAXMU=0.05
MINBETA=0.0
MAXBETA=0.9
MINPGEOM=0.7
MAXPGEOM=1.0
SCALE=0.4
GAMMA=1.2