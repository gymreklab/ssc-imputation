startchrom=1
endchrom=22
VCFDIR=/oasis/projects/nsf/csd568/mgymrek/ssc-quads/hipstr_vcfs/final
OUTDIR=/oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/all_mutations
FINALOUTDIR=/oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered_mutations
LOGDIR=/oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered_mutations/logs
FILTERDIR=/oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered_mutations/filter_locus
ANNDIR=/oasis/projects/nsf/csd568/ileena/ssc_phase1_denovos/filtered_mutations/annotations

# Params for STRDenovoTools                                                                                                            
PTHRESH=0.8
MAXALLELES=100
MINCOV=10
MINSPANCOV=10
MINSCORE=0.9
MINSUPPREADS=2
MAXCOUNT=10 # for rare analysis                                                                                                        

# Annotation and filtering files                                                                                                       
CODING=/oasis/projects/nsf/csd568/mgymrek/data/coding.bed
UTR3=/oasis/projects/nsf/csd568/mgymrek/data/3utr.bed
UTR5=/oasis/projects/nsf/csd568/mgymrek/data/5utr.bed
INTRON=/oasis/projects/nsf/csd568/mgymrek/data/introns.bed
DONOR=/oasis/projects/nsf/csd568/mgymrek/data/donorsites.bed
ACCEPTOR=/oasis/projects/nsf/csd568/mgymrek/data/acceptorsites.bed
PROMOTER1KB=/oasis/projects/nsf/csd568/mgymrek/data/hg19_promoter_1kb.bed
PROMOTER3KB=/oasis/projects/nsf/csd568/mgymrek/data/hg19_promoter_3kb.bed
PROMOTER5KB=/oasis/projects/nsf/csd568/mgymrek/data/hg19_promoter_5kb.bed
EXAC=/oasis/projects/nsf/csd568/mgymrek/data/fordist_cleaned_exac_nonTCGA_z_pli_rec_null_data.txt
GENES=/oasis/projects/nsf/csd568/mgymrek/data/hg19_refseq_genes.bed
CONSTRAINT=/oasis/projects/nsf/csd568/mgymrek/data/ssc_autosomal_perlocus_constraint.bed
HIPPROP=/oasis/projects/nsf/csd568/mgymrek/data/GRCh37.hipstr_reference_sorted_properties.tab.gz
ASDGENES=/oasis/projects/nsf/csd568/mgymrek/data/SFARI-Gene_genes_export13-11-2017.csv
RNABP=/oasis/projects/nsf/csd568/mgymrek/data/hipstr_rnabp_hg19_total.bed
RNABP2=/oasis/projects/nsf/csd568/mgymrek/data/hipstr_rnabp_hg19_TARDBP.bed
TF=/oasis/projects/nsf/csd568/mgymrek/data/hipstr_tfbs_hg19_total.bed
TFALL=/oasis/projects/nsf/csd568/mgymrek/data/hipstr_tfbs_hg19.bed
HISTONE=/oasis/projects/nsf/csd568/mgymrek/data/hipstr_histone_hg19.bed
ESTRS=/oasis/projects/nsf/csd568/mgymrek/data/GTEx_merged_causality.tab
CYTOBANDS=/oasis/projects/nsf/csd568/mgymrek/data/cytoBand.txt.gz
ASDCNV=/oasis/projects/nsf/csd568/mgymrek/data/SFARI-CNV-annotated.tab
BRAIN=/oasis/projects/nsf/csd568/mgymrek/data/Brain_expressed_tissues.bed

# Annotation and filtering parameters                                                                                                  
MINCHILDREN=50
MAXFILTFAM=5
