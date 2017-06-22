#!/usr/bin/env python
"""
Calculate significance of de novo mutation excess

Usage:
./calc_denovo_sig.py <mutation counts>
"""

import pandas as pd
from scipy import stats
import numpy as np
import sys

try:
    datafile = sys.argv[1]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

data = pd.read_csv(datafile, sep="\t", names=["chrom","pos","total", \
                                              "mut_affected_length", "mut_unaffected_length", \
                                              "mut_affected_seq","mut_unaffected_seq"], usecols=range(7))

# overall - length
total = sum(data["total"])
table_length = np.array([[sum(data["mut_affected_length"]), total-sum(data["mut_affected_length"])], \
                  [sum(data["mut_unaffected_length"]), total-sum(data["mut_unaffected_length"])]])
res_length = list(stats.fisher_exact(table_length, alternative="greater"))
sys.stdout.write("\t".join(map(str, ["# Overall - length"] + res_length))+"\n")

# overall - seq
table_seq = np.array([[sum(data["mut_affected_seq"]), total-sum(data["mut_affected_seq"])], \
                  [sum(data["mut_unaffected_seq"]), total-sum(data["mut_unaffected_seq"])]])
res_seq = list(stats.fisher_exact(table_seq, alternative="greater"))
sys.stdout.write("\t".join(map(str, ["# Overall - seq"] + res_seq))+"\n")

# per locus
def GetDenovoExcess(x):
    table = np.array([[x["mut_affected_length"], x["total"]-x["mut_affected_length"]], \
                      [x["mut_unaffected_length"], x["total"]-x["mut_unaffected_length"]]])
    return stats.fisher_exact(table, alternative="greater")[1]
def GetDenovoExcessSeq(x):
    table = np.array([[x["mut_affected_seq"], x["total"]-x["mut_affected_seq"]], \
                      [x["mut_unaffected_seq"], x["total"]-x["mut_unaffected_seq"]]])
    return stats.fisher_exact(table, alternative="greater")[1]
                     
data["pval_length"] = data.apply(GetDenovoExcess, 1)
data["pval_seq"] = data.apply(GetDenovoExcessSeq, 1)
data[["chrom","pos","pval_length","pval_seq", \
      "total","mut_affected_length","mut_unaffected_length", \
      "mut_affected_seq","mut_unaffected_seq"]].sort_values("pval_length").to_csv(sys.stdout, sep="\t", index=False)
