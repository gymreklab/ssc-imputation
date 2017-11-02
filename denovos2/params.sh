startchrom=1
endchrom=22
OUTDIR=/storage/mgymrek/ssc-denovos/denovos2/denovocalls
FINALOUTDIR=/storage/mgymrek/ssc-denovos/denovos2/denovocalls_filtered
LOGDIR=/storage/mgymrek/ssc-denovos/denovos2/denovocalls/log
FILTERDIR=/storage/mgymrek/ssc-denovos/denovos2/denovocalls/filters
ANNDIR=/storage/mgymrek/ssc-denovos/denovos2/denovocalls/annotations

# Params for STRDenovoTools
PTHRESH=0.8
MAXALLELES=100
MINCOV=10
MINSPANCOV=10
MINSCORE=0.9
MINSUPPREADS=2

# Annotation and filtering files
SEGDUP=/storage/resources/dbase/human/hg19/hg19_segmentalduplications.bed
CODING=/storage/resources/dbase/human/hg19/hg19_codingexons.bed
EXAC=/storage/mgymrek/ssc-denovos/denovos2/other-data/fordist_cleaned_exac_nonTCGA_z_pli_rec_null_data.txt
GENES=/storage/mgymrek/ssc-denovos/denovos2/other-data/hg19_refseq_genes.bed
CONSTRAINT=/storage/mgymrek/ssc-denovos/mutea-results/ssc_autosomal_perlocus_constraint.bed
HIPPROP=/storage/mgymrek/ssc-denovos/denovos2/other-data/GRCh37.hipstr_reference_sorted_properties.tab.gz

# Annotation and filtering parameters
MINCHILDREN=50
MAXFILTFAM=5
