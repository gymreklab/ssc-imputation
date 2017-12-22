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
MAXCOUNT=10 # for rare analysis

# Annotation and filtering files
SEGDUP=/storage/resources/dbase/human/hg19/hg19_segmentalduplications.bed
CODING=/storage/mgymrek/gtex/annotations/coding.bed
UTR3=/storage/mgymrek/gtex/annotations/3utr.bed
UTR5=/storage/mgymrek/gtex/annotations/5utr.bed
INTRON=/storage/mgymrek/gtex/annotations/introns.bed
DONOR=/storage/mgymrek/gtex/annotations/donorsites.bed
ACCEPTOR=/storage/mgymrek/gtex/annotations/acceptorsites.bed
PROMOTER1KB=/storage/mgymrek/gtex/annotations/hg19_promoter_1kb.bed
PROMOTER3KB=/storage/mgymrek/gtex/annotations/hg19_promoter_3kb.bed
PROMOTER5KB=/storage/mgymrek/gtex/annotations/hg19_promoter_5kb.bed
EXAC=/storage/mgymrek/ssc-denovos/denovos2/other-data/fordist_cleaned_exac_nonTCGA_z_pli_rec_null_data.txt
GENES=/storage/mgymrek/ssc-denovos/denovos2/other-data/hg19_refseq_genes.bed
CONSTRAINT=/storage/mgymrek/ssc-denovos/mutea-results/ssc_autosomal_perlocus_constraint.bed
HIPPROP=/storage/mgymrek/ssc-denovos/denovos2/other-data/GRCh37.hipstr_reference_sorted_properties.tab.gz
ASDGENES=/storage/mgymrek/ssc-denovos/denovos2/other-data/SFARI-Gene_genes_export13-11-2017.csv
RNABP=/storage/mgymrek/gtex/causality/features/hipstr_rnabp_hg19_total.bed
RNABP2=/storage/mgymrek/gtex/causality/features/hipstr_rnabp_hg19_TARDBP.bed
TF=/storage/mgymrek/gtex/causality/features/hipstr_tfbs_hg19_total.bed
TFALL=/storage/mgymrek/gtex/causality/features/hipstr_tfbs_hg19.bed
HISTONE=/storage/mgymrek/gtex/causality/features/hipstr_histone_hg19.bed
ESTRS=/storage/mgymrek/gtex/causality/GTEx_merged_causality.tab 
CYTOBANDS=/storage/mgymrek/ssc-denovos/denovos2/other-data/cytoBand.txt.gz
ASDCNV=/storage/mgymrek/ssc-denovos/denovos2/other-data/SFARI-CNV-annotated.tab
BRAIN=/storage/mgymrek/ssc-denovos/denovos2/other-data/Brain_expressed_tissues.bed 

# Annotation and filtering parameters
MINCHILDREN=50
MAXFILTFAM=5
