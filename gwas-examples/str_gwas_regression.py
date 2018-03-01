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

def RunRegression(data):
    formula = "phenotype ~ gtsum+C(sex)+C(cohort)+"+"+".join(COVARCOLS)
    pgclogit = logit(formula=formula, data=data).fit(disp=0)
    return pgclogit

def PrintLine(chrom, pos, testclass, res, nsamp, f):
    f.write("\t".join(map(str, [chrom, pos, testclass, res.params["gtsum"], res.pvalues["gtsum"], nsamp]))+"\n")
    f.flush()

def GetAlleles(gts):
    alleles = set()
    for item in gts.values:
        alleles.add(item[0])
        alleles.add(item[1])
    return list(alleles)

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--strvcf", help="VCF file with STR genotypes", required=True, type=str)
    parser.add_argument("--sampledata", help="Tab file with sample,phenotype for each sample", required=True, type=str)
    parser.add_argument("--out", help="Prefix for output files. Default: stdout", required=False, type=str)
    parser.add_argument("--region", help="Restrict analysis to STRs in this region.", type=str)
    args = parser.parse_args()

    # Prepare output
    if args.out: outf = args.out + ".gtsum.tab"
    else: outf = sys.stdout

    # Load sample metadata
    sdata = pd.read_csv(args.sampledata, sep="\t", names=["sample","phenotype"])
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
        data["phenotype"] = data["phenotype"]-1
        # Run on sum of allele lengths
        data["gtsum"] = data["gts"].apply(lambda x: sum(x))
        try:
            res = RunRegression(data)
        except: continue
        PrintLine(record.CHROM, record.POS, "locus", res, data.shape[0], outf)
        # Run on each allele separately
#        alleles = GetAlleles(data["gts"])
#        for al in alleles:
#            data["gtsum"] = data.apply(lambda x: int(x["gts"][0]==al) + int(x["gts"][1]==al), 1)
#            try:
#                res = RunRegression(data)
#            except: continue
#            PrintLine(record.CHROM, record.POS, al, res, data.shape[0], outf)
        
if __name__ == "__main__":
    main()
