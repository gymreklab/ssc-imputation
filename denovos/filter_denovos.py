#!/usr/bin/env python
"""
Filter processed de novo calls
"""

import argparse
import numpy as np
import pandas as pd
import scipy.misc
import sys

# Genotype info columns
q_col = 1
dstutter_col = 2
dp_col = 3
dflankindel_col = 4
mallreads_col = 5
gb_col = 6

# TODO check these calculations
def GetScore(nomut_ll, anymut_ll, \
             denovo_ll_affected, denovo_ll_unaffected, \
             other_ll_affected, other_ll_unaffected):
    """
    P(i|D) = P(D|i)P(i)/[sum_i P(D|i)P(i)]
    P(D|i) are ll scores

    For nomut vs. anymut
    P(anymut) = mutation rate
    P(mut) = 1-mutation rate

    For nomut, mut1, mut2
    mut1 = denovo1 + other1
    mut2 = denovo2 + other2
    p(nomut) = 1-2*mutrate
    P(mut1) = p(mut2) = mutrate
    return posterior prob of nomut, mut1, mut2
    """
    # Combine "other" and "denovo" probs for each child
    # Convert to natural logs
    base = np.log10(np.exp(1))
    mut_affected_ll = scipy.misc.logsumexp([denovo_ll_affected/base, other_ll_affected/base])
    mut_unaffected_ll = scipy.misc.logsumexp([denovo_ll_unaffected/base, other_ll_unaffected/base])
    nomut_ll = nomut_ll/base

    # Get denom
    p_mut = 0.0001 # TODO use per locus
    denom_elts = [nomut_ll+np.log(1-2*p_mut), \
                  mut_affected_ll+np.log(p_mut), \
                  mut_unaffected_ll+np.log(p_mut)]
    denom = scipy.misc.logsumexp(denom_elts)

    # Get posteriors
    post_nomut = np.exp(nomut_ll+np.log(1-2*p_mut) - denom)
    post_mut_affected = np.exp(mut_affected_ll+np.log(p_mut) - denom)
    post_mut_unaffected = np.exp(mut_unaffected_ll+np.log(p_mut) - denom)
    return post_nomut, post_mut_affected, post_mut_unaffected

def GetSpanning(mallreads):
    items = mallreads.split(";")
    total = 0
    for item in items:
        total += int(item.split("|")[1])
    return total

def GetMinPercSupp(mallreads, gb):
    try:
        alleles = map(int, gb.split("|"))
    except: alleles = map(int, gb.split("/"))
    mallreads_dict = {}
    for item in mallreads.split(";"):
        allele, dp = map(int, item.split("|"))
        mallreads_dict[allele] = mallreads_dict.get(allele, 0) + dp
    total = sum(mallreads_dict.values())
    perc_allele1 = mallreads_dict.get(alleles[0], 0)*1.0/total
    perc_allele2 = mallreads_dict.get(alleles[1], 0)*1.0/total
    return min([perc_allele1, perc_allele2])

def IsUnit(gb, period):
    alleles = map(int, gb.split("|"))
    if alleles[0]%period != 0: return False
    if alleles[1]%period != 0: return False
    return True

def GetFamilyItems(x):
    return {
        "affected": x["gtdata_affected"].split(":"),
        "unaffected": x["gtdata_unaffected"].split(":"),
        "father": x["gtdata_father"].split(":"),
        "mother": x["gtdata_mother"].split(":")
    }

def FilterIndividual(gtitems, minq, minst, min_span, min_percsupp, unit, period):
    if gtitems[mallreads_col] == ".": return True
    if float(gtitems[q_col]) < minq: return True
    if float(gtitems[dstutter_col])*1.0/float(gtitems[dp_col]) > minst: return True
    if float(gtitems[dflankindel_col])*1.0/float(gtitems[dp_col]) > minst: return True
    total_span = GetSpanning(gtitems[mallreads_col])
    if total_span < min_span: return True
    percsupp = GetMinPercSupp(gtitems[mallreads_col], gtitems[gb_col])
    if percsupp < min_percsupp: return True
    if unit and not IsUnit(gtitems[gb_col], period): return True
    return False

def FilterFamily(x, minq_child, minq_father, minq_mother, \
                   minst_child, minst_father, minst_mother, \
                   min_span_child, min_span_father, min_span_mother, \
                   min_percsupp_child, min_percsupp_father, min_percsupp_mother, unit):
    # Get items
    family_items = GetFamilyItems(x)
    # Individual level filters
    if FilterIndividual(family_items["affected"], \
                        minq_child, minst_child, \
                        min_span_child, min_percsupp_child, \
                        unit, x["period"]): return True
    if FilterIndividual(family_items["unaffected"], \
                        minq_child, minst_child, \
                        min_span_child, min_percsupp_child, \
                        unit, x["period"]): return True
    if FilterIndividual(family_items["mother"], \
                        minq_mother, minst_mother, \
                        min_span_mother, min_percsupp_mother, \
                        unit, x["period"]): return True
    if FilterIndividual(family_items["father"], \
                        minq_father, minst_father, \
                        min_span_father, min_percsupp_father, \
                        unit, x["period"]): return True
    return False

def GetFamilyGB(x):
    family_items = GetFamilyItems(x)
    gb = {}
    for key in family_items:
        gb[key] = map(int, family_items[key][gb_col].split("|"))
    return gb

def GetMendLengths(x):
    """
    0: neither sibling follows mendelian
    1: one does one doesn't
    2: both siblings follow mendelian inheritance
    """
    family_gb = GetFamilyGB(x)
    # Test each child
    affected = False
    unaffected = False
    if family_gb["affected"][0] in family_gb["father"] and \
       family_gb["affected"][1] in family_gb["mother"]: affected = True
    if family_gb["affected"][1] in family_gb["father"] and \
       family_gb["affected"][0] in family_gb["mother"]: affected = True
    if family_gb["unaffected"][0] in family_gb["father"] and \
       family_gb["unaffected"][1] in family_gb["mother"]: unaffected = True
    if family_gb["unaffected"][1] in family_gb["father"] and \
       family_gb["unaffected"][0] in family_gb["mother"]: unaffected = True
    return int(affected)+int(unaffected)

def GetMendLengthsOther(x):
    """
    If lengths of genotypes can be explained by mend inheritance with dropout
    """
    family_gb = GetFamilyGB(x)
    parent_alleles = family_gb["mother"] + family_gb["father"]
    affected = False
    unaffected = False
    if family_gb["affected"][0] in parent_alleles and \
       family_gb["affected"][1] in parent_alleles: affected = True
    if family_gb["unaffected"][0] in parent_alleles and \
       family_gb["unaffected"][1] in parent_alleles: unaffected = True
    return int(affected)+int(unaffected)

def GetNewAllele(x):
    """
    Return new allele in child genotypes. Assume only one
    """
    family_gb = GetFamilyGB(x)
    parent_alleles = family_gb["mother"] + family_gb["father"]
    child_alleles = family_gb["affected"] + family_gb["unaffected"]
    for al in child_alleles:
        if al not in parent_alleles: return al
    return "NA"

def GetSupportingReads(mallreads, testallele):
    for item in mallreads.split(";"):
        allele, dp = map(int, item.split("|"))
        if allele == testallele: return dp
    return 0

def GetDeNovoFilter(x, parent_allele_count):
    """
    For potential de novo calls, return True if:
    new allele found at >parent_allele_count reads in one of parents
    """
    if x["parent_origin"] == "NA": return False
    if x["filter"]: return False
    family_items = GetFamilyItems(x)
    supp_mother = GetSupportingReads(family_items["mother"][mallreads_col], x["new_allele"])
    supp_father = GetSupportingReads(family_items["father"][mallreads_col], x["new_allele"])
    return (supp_mother > parent_allele_count) or \
        (supp_father > parent_allele_count)

def GetPOO(x):
    """
    For de novo allele, return mother or father or both
    If no de novo allele, return NA
    """
    if x["mend_lengths"] or x["mend_lengths_other"] or x["new_allele"] == "NA": return "NA"
    new_allele = x["new_allele"]
    family_gb = GetFamilyGB(x)
    poo = None
    if new_allele in family_gb["affected"]:
        old_allele = [item for item in family_gb["affected"] if item != new_allele]
        try:
            old_allele = old_allele[0]
            if old_allele in family_gb["mother"] and old_allele in family_gb["father"]:
                return "both"
            if old_allele in family_gb["mother"] and old_allele not in family_gb["father"]:
                return "father"
            if old_allele in family_gb["father"] and old_allele not in family_gb["mother"]:
                return "mother"
        except: pass
    if new_allele in family_gb["unaffected"]:
        old_allele = [item for item in family_gb["unaffected"] if item != new_allele]
        try:
            old_allele = old_allele[0]
            if old_allele in family_gb["mother"] and old_allele in family_gb["father"]:
                return "both"
            if old_allele in family_gb["mother"] and old_allele not in family_gb["father"]:
                return "father"
            if old_allele in family_gb["father"] and old_allele not in family_gb["mother"]:
                return "mother"
        except: pass
    return "NA"

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--infile", help="Input file", type=str, required=True)
    parser.add_argument("--outfile", help="Output file", type=str, required=True)
    parser.add_argument("--minq-child", help="Min GT qual score for child", type=float, default=0.9)
    parser.add_argument("--minq-father", help="Min GT qual score for father", type=float, default=0.9)
    parser.add_argument("--minq-mother", help="Min GT qual score for mother", type=float, default=0.9)
    parser.add_argument("--minst-child", help="Min % stutter/indel reads for child", type=float, default=0.2)
    parser.add_argument("--minst-father", help="Min % stutter/indel reads for father", type=float, default=0.2)
    parser.add_argument("--minst-mother", help="Min % stutter/indel reads for mother", type=float, default=0.2)
    parser.add_argument("--min-spanning-child", help="Min number of spanning reads for child", type=int, default=10)
    parser.add_argument("--min-spanning-father", help="Min number of spanning reads for father", type=int, default=10)
    parser.add_argument("--min-spanning-mother", help="Min number of spanning reads for mother", type=int, default=10)
    parser.add_argument("--min-percsupp-child", help="Min number of reads supporting alleles for child", type=float, default=0.2)
    parser.add_argument("--min-percsupp-father", help="Min number of reads supporting alleles  for father", type=float, default=0.2)
    parser.add_argument("--min-percsupp-mother", help="Min number of reads supporting alleles for mother", type=float, default=0.2)
    parser.add_argument("--parent-allele-count", help="Set denovo filter if > this many reads support new allele in a parent", type=int, default=0)
    parser.add_argument("--unit", help="Require bp changes to be unit number of repeats", action="store_true")
    args = parser.parse_args()

    outcols = ["chrom","pos","period", \
               "post_nomut","post_mut_affected","post_mut_unaffected", \
               "mend_lengths", "mend_lengths_other", \
               "gtdata_affected","gtdata_unaffected","gtdata_mother","gtdata_father", \
               "filter","denovo_filter","parent_origin","new_allele","ref","alt"]

    # Load input file
    data = pd.read_csv(args.infile, sep="\t")

    # Remove lines with missing denovo data
    data = data[data["nomut_ll"].apply(str) != "."]
    if data.shape[0] == 0:
        f = open(args.outfile, "w")
        f.write("\t".join(outcols)+"\n")
        f.close()
        sys.exit(1)

    # Make de novo score
    data["nomut_ll"] = data["nomut_ll"].apply(float)
    data["anymut_ll"] = data["anymut_ll"].apply(float)
    data["denovo_ll_affected"] = data["denovo_ll"].apply(lambda x: float(x.split(",")[0]))
    data["denovo_ll_unaffected"] = data["denovo_ll"].apply(lambda x: float(x.split(",")[1]))
    data["other_ll_affected"] = data["other_ll"].apply(lambda x: float(x.split(",")[0]))
    data["other_ll_unaffected"] = data["other_ll"].apply(lambda x: float(x.split(",")[1]))
    scores = data.apply(lambda x: GetScore(x["nomut_ll"], x["anymut_ll"], \
                                           x["denovo_ll_affected"], x["denovo_ll_unaffected"], \
                                           x["other_ll_affected"], x["other_ll_unaffected"]), 1)
    data["post_nomut"] = scores.apply(lambda x: x[0])
    data["post_mut_affected"] = scores.apply(lambda x: x[1])
    data["post_mut_unaffected"] = scores.apply(lambda x: x[2])

    # Implement filters
    data["filter"] = data.apply(lambda x: \
                                    FilterFamily(x, args.minq_child, args.minq_father, args.minq_mother, \
                                                   args.minst_child, args.minst_father, args.minst_mother, \
                                                   args.min_spanning_child, args.min_spanning_father, args.min_spanning_mother, \
                                                   args.min_percsupp_child, args.min_percsupp_father, args.min_percsupp_mother, \
                                                   args.unit), 1)

    data["mend_lengths"] = data.apply(lambda x: GetMendLengths(x), 1)
    data["mend_lengths_other"] = data.apply(lambda x: GetMendLengthsOther(x), 1)
    data["new_allele"] = data.apply(lambda x: GetNewAllele(x), 1)
    data["parent_origin"] = data.apply(lambda x: GetPOO(x), 1)
    data["denovo_filter"] = data.apply(lambda x: GetDeNovoFilter(x, args.parent_allele_count), 1)
    data[outcols].sort_values("post_nomut").to_csv(args.outfile, sep="\t", index=False)

if __name__ == "__main__":
    main()
