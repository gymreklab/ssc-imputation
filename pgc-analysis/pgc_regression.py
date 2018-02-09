#!/usr/bin/env python
"""
Regression analysis between STRs and SCZ phenotype in PGC
"""

import argparse
import numpy as np
import pandas as pd
from statsmodels.formula.api import logit
import sys
import vcf

COVARCOLS = ["C1","C2","C3","C4","C5","C6","C7","C9","C15","C18"]

def RunRegression(data):
    formula = "phenotype ~ gtsum+"+"+".join(COVARCOLS+["sex"])
    pgclogit = logit(formula=formula, data = data[["phenotype","gtsum"]+COVARCOLS+["sex"]]).fit(disp=0)
    return pgclogit

def PrintLine(chrom, pos, res, f):
    f.write("\t".join(map(str, [chrom, pos, res.params["gtsum"], res.pvalues["gtsum"]]))+"\n")

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--strvcf", help="VCF file with STR genotypes", required=True, type=str)
    parser.add_argument("--sampledata", help="Tab file with pheno/covars for each sample", required=True, type=str)
    parser.add_argument("--out", help="Prefix for output files. Default: stdout", required=False, type=str)
    parser.add_argument("--region", help="Restrict analysis to STRs in this region.", type=str)
    args = parser.parse_args()

    # Prepare output
    if args.out: outf = args.out + ".gtsum.tab"
    else: outf = sys.stdout

    # Load sample metadata
    sdata = pd.read_csv(args.sampledata, sep="\t", names=["sample","sex","phenotype","cohort"] + COVARCOLS)
    str_reader = vcf.Reader(open(args.strvcf, "rb"))
    str_records = str_reader
    if args.region:
        chrom = args.region.split(":")[0]
        start = int(args.region.split(":")[1].split("-")[0])
        end = int(args.region.split(":")[1].split("-")[1])
        str_records = str_reader.fetch(chrom, start, end)

    for record in str_records:
        samples = []
        gts = []
        if None in record.ALT: continue # not polymorphic
        allelelens = [0] + [len(item)-len(record.REF) for item in record.ALT]
        for sample in record:
            samples.append(sample.sample)
            gts.append(map(lambda x: allelelens[int(x)], sample.gt_alleles))
        gdata = pd.DataFrame({"sample": samples, "gts": gts})
        data = pd.merge(gdata, sdata, on=["sample"])
        data["gtsum"] = data["gts"].apply(lambda x: sum(x))
        data["phenotype"] = data["phenotype"]-1
        res = RunRegression(data)
        PrintLine(record.CHROM, record.POS, res, outf)

if __name__ == "__main__":
    main()
