# Get SSC stats
./preprocess_wrapper.sh
./combine_stats.sh

# Get lobSTR catalog stats
~/workspace/mgymrek-utils/vcf_het.py \
    --vcf /storage/mgymrek/ssc-imputation/tmp/phase_1_final_calls.vcf.gz > \
    /storage/mgymrek/ssc-imputation/tmp/lobstr_het_stats.tab
