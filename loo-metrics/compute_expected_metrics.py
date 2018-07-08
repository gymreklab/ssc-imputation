#!/usr/bin/env python

"""
Compute LD metrics (allele r2, lengthr2)
for each STR as expected at random.

For each STR, for "true" genotypes, randomly generate
genotypes based on observed allele frequencies
"""

import argparse
import numpy as np
import random
import scipy.stats
import sys
import vcf

def GetGTs(record):
    gts = []
    for sample in record:
        if sample["GB"] is not None:
            a1, a2 = sample["GB"].split("|")
            gts.append((int(a1), int(a2)))
    return gts

def PermuteGTs(obs_gts):
    alleles = []
    for a1, a2 in obs_gts:
        alleles.extend([a1, a2])
    random.shuffle(alleles)
    gts = []
    for i in range(len(alleles)/2):
        gts.append((alleles[i*2], alleles[i*2+1]))
    return gts

def GetLengthR2(gts1, gts2):
    alleles1 = [item[0] for item in gts1]+[item[1] for item in gts1]
    alleles2 = [item[0] for item in gts2]+[item[1] for item in gts2]
    return scipy.stats.pearsonr(alleles1, alleles2)[1]**2

def GetAlleleR2(gts1, gts2):
    alleles = set([item[0] for item in gts1]+[item[1] for item in gts2])
    res = []
    for a in alleles:
        gts1_a = [(int(item[0]==a), int(item[1]==a)) for item in gts1]
        gts2_a = [(int(item[0]==a), int(item[1]==a)) for item in gts2]
        gts1_a = [item[0] for item in gts1_a]+[item[1] for item in gts1_a]
        gts2_a = [item[0] for item in gts2_a]+[item[1] for item in gts2_a]
        a_r2 = scipy.stats.pearsonr(gts1_a, gts2_a)[1]**2
        res.append((a, a_r2))
    return res

def GetConc(gts1, gts2):
    return np.mean([(gts1[i][0]==gts2[i][0] and gts1[i][1]==gts2[i][1]) for i in range(len(gts1))])

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--vcf", help="VCF file with true genotypes", type=str, required=True)
    parser.add_argument("--out", help="Output file prefix", type=str, required=True)
    args = parser.parse_args()

    f_locus = open(args.out + ".locus.tab", "w")
    f_locus.write("\t".join(["chrom","pos","start","length_r2","conc"])+"\n")
    f_locus.flush()
    f_allele = open(args.out + ".allele.tab", "w")
    f_allele.write("\t".join(["chrom","pos","start","allele","allele_r2"])+"\n")
    f_allele.flush()
    reader = vcf.Reader(open(args.vcf, "rb"))
    for record in reader:
        if len(record.FILTER)!=0: continue
        obs_gts = GetGTs(record)
        perm_gts = PermuteGTs(obs_gts)
        length_r2 = GetLengthR2(obs_gts, perm_gts)
        conc = GetConc(obs_gts, perm_gts)
        # Locus metrics
        f_locus.write("\t".join(map(str, [record.CHROM, record.POS, record.INFO["START"], \
                                             length_r2, conc]))+"\n")
        f_locus.flush()
        # Allele metrics
        allele_r2s = GetAlleleR2(obs_gts, perm_gts)
        sys.exit(1)
        for allele, allele_r2 in allele_r2s:
            f_allele.write("\t".join(map(str, [record.CHROM, record.POS, record.INFO["START"], \
                                                 allele, allele_r2]))+"\n")
            f_allele.flush()
if __name__ == "__main__":
    main()
