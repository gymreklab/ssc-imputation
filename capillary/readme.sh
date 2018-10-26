#!/bin/bash


# Get list of loci to use
cat marshfield_marker_corrections.tab | sort -k 4,4 > marshfield_marker_corrections_sorted.tab
cat /storage/s1saini/manuscript_strsnp/fig3/cap_elec/Pemberton_AdditionalFile1_11242009.txt | cut -f 4,5 | sort -k 5,5 > pemberton_sorted.tab
