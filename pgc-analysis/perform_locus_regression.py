#!/usr/bin/env python2.7

"""
Perform PGC Regression. Logistic regression model:

phenotype ~ genotype + cohort + sex + PCs TODO change these?

Where genotype is either:
1. sum of STR allele lengths
2. Separate model for each allele treating it as biallelic
3. Each SNP position

Usage:
./perform_atxn7_regression.py --data <regression data.tab>
"""

import argparse
import numpy as np
import pandas as pd
from statsmodels.formula.api import logit
import sys

COVARCOLS = ["C1", "C2", "C3", "C4", "C5", "C6", "C7", "C9", "C15", "C18", "cohort"] # add back sex? 
MINMAF=0.001

def RunAndPrint(data, name, maf):
    formula = "phenotype ~ genotype+"+"+".join(COVARCOLS)
    try:
        pgclogit = logit(formula=formula, data=data[["phenotype", "genotype"]+COVARCOLS]).fit(disp=0)
    except:
        sys.stderr.write("ERROR running regression for %s, maf=%s\n"%(name, maf))
        return
    coef = pgclogit.params["genotype"]
    pval = pgclogit.pvalues["genotype"]
    sys.stdout.write("\t".join(map(str, [name, coef, pval, maf]))+"\n")
    sys.stdout.flush()

def GetMAF(genotypes):
    total = len(genotypes)*2
    return sum(genotypes)*1.0/total

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--data", help="Input regression data, tab delim", type=str, required=True)
    args = parser.parse_args()

    # Output
    sys.stdout.write("\t".join(["Position","coefficient","pval","maf"])+"\n")
    sys.stdout.flush()
    # Clean up the data frame
    data = pd.read_csv(args.data, sep="\t")
    data["phenotype"] = data["phenotype"]-1

    # Regress on sum of alleles
    data["genotype"] = data["strgt1"] + data["strgt2"]
    RunAndPrint(data, "STR", "-1")

    # Regress on max allele
    data["genotype"] = data.apply(lambda x: max([x["strgt1"], x["strgt2"]]), 1)
    RunAndPrint(data, "STR-max", "-1")

    # Regress on allele diff
    data["genotype"] = data.apply(lambda x: abs(x["strgt1"]-x["strgt2"]), 1)
    RunAndPrint(data, "STR-diff", "-1")

    # Regress on each allele separately
    alleles = set(list(data["strgt1"].values)+list(data["strgt2"].values))
    for al in alleles:
        data["genotype"] = data.apply(lambda x: int(x["strgt1"]==al)+int(x["strgt2"]==al), 1)
        RunAndPrint(data, "STR-%s"%al, GetMAF(list(data["genotype"])))

    # Regress on each SNP
    snpcols = [item for item in list(data.columns) if item.startswith("SNP_")]
    for item in snpcols:
        data["genotype"] = data[item]
        maf = GetMAF(data["genotype"])
        if maf < MINMAF or (1-maf) < MINMAF: continue
        RunAndPrint(data, item, maf)

if __name__ == "__main__":
    main()
