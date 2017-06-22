# TODO change ./test.sh runs denovofinder. change to run on all hipstr results

./process_families.sh puts denovofinder res per family into tab format with fields for filtering. Calls:
  - ./process_denovo.sh mostly wrapper to call filter_denovos.py script
  - ./filter_denovos.py applies filters
  - ./summarize_family.py summarizes family info

./combine_family_results.sh compiles all family results into a single file
  - calc_denovo_sig.py calculates sig of excess in affected vs. unaffected perlocus, and overall
