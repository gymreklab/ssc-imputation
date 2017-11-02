#!/bin/bash

source params.sh

tmpdir=$(mktemp -d)

# Motif, length, uninterrupted track length
./get_motif_length_features.py ${HIPREF} ${REFFA} > ${tmpdir}/hipstr_ref_features.bed

# Combine (might add more features later)
echo "chrom,start,end,motif,length,uninterrupted_length" | \
    sed 's/,/\t/g' > ${OUTDIR}/GRCh37.hipstr_reference_sorted_properties.tab
cat ${tmpdir}/hipstr_ref_features.bed >> ${OUTDIR}/GRCh37.hipstr_reference_sorted_properties.tab
gzip -f ${OUTDIR}/GRCh37.hipstr_reference_sorted_properties.tab
