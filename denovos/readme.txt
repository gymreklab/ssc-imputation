./run_denovofinder.sh runs DenovoFinder on all chromosomes
./get_locus_stats.sh gets stats for each locus

./process_families.sh puts denovofinder res per family into tab format with fields for filtering. Calls:
  - ./process_denovo.sh mostly wrapper to get ready to call filter_denovos.py
  - ./filter_denovos.py applies filters
  - ./summarize_family.py summarizes family info

./combine_family_results.sh compiles all family results into a single file
  - calc_denovo_sig.py calculates sig of excess in affected vs. unaffected perlocus, and overall
