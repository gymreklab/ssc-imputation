#!/usr/bin/env python
"""
Using filtered denovo calls, output stats per locus
"""

import argparse
import os
import pandas as pd
import scipy.stats
import sys
import tabix

# For profiling
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput, GephiOutput

def MSG(msg):
    sys.stderr.write(msg.strip() + "\n")

def ERROR(msg):
    MSG(msg)
    sys.exit(1)

def ProcessLocus(df, pthresh, outm, outmcols, filter_both_kids, max_filtered_families, min_children):
    # Apply filters first
    if filter_both_kids:
        mutfam = df[(df["posterior"]>=pthresh)].groupby("family", as_index=False).agg({"child": len})
        filterfam = list(mutfam[mutfam["child"]>1]["family"].values)
        if len(filterfam) > max_filtered_families:
            MSG("Families %s for locus %s:%s fail filter. Skipping locus"%(",".join(map(str, filterfam)), df["chrom"].values[0], df["pos"].values[0]))
            return None
        if len(filterfam) > 0:
            MSG("Removing families %s for locus %s:%s"%(",".join(map(str, filterfam)), df["chrom"].values[0], df["pos"].values[0]))
            df = df[df["family"].apply(lambda x: x not in filterfam)]
    if df.shape[0] == 0:
        MSG("No children found")
        return None
    if df[(df["phenotype"]==1)].shape[0] < min_children:
        MSG("Not enough children for locus %s:%s - %s"%(df["chrom"].values[0], df["pos"].values[0], df.shape[0]))
        return None
    period = df["period"].values[0]
    total_children = df.shape[0]
    total_mutations = df[(df["posterior"]>=pthresh)].shape[0]
    total_mutation_rate = total_mutations*1.0/total_children
    affected_children = df[(df["phenotype"]==2)].shape[0]
    affected_mutations = df[(df["phenotype"]==2) & (df["posterior"]>=pthresh)].shape[0]
    affected_mutation_rate = affected_mutations*1.0/affected_children
    unaffected_children = df[(df["phenotype"]==1)].shape[0]
    unaffected_mutations = df[(df["phenotype"]==1) & (df["posterior"]>=pthresh)].shape[0]
    unaffected_mutation_rate = unaffected_mutations*1.0/unaffected_children
    n11 = affected_mutations
    n12 = affected_children - affected_mutations
    n21 = unaffected_mutations
    n22 = unaffected_children - unaffected_mutations
    pvalue = scipy.stats.fisher_exact([[n11,n12],[n21, n22]], alternative="greater")[1]
    children_with_mutations = ",".join(list(df[(df["posterior"]>=pthresh)].apply(lambda x: x["family"]+":"+x["child"], 1).values))
    # Output mutations
    df[(df["posterior"]>=pthresh)][["chrom"]+outmcols[1:]].to_csv(outm, header=False, index=False, sep="\t")
    outm.flush()

    # Return summary line
    return [period, \
        total_children, total_mutations, total_mutation_rate, \
        affected_children, affected_mutations, affected_mutation_rate, \
        unaffected_children, unaffected_mutations, unaffected_mutation_rate, \
        pvalue, children_with_mutations]
        

def ReadDenovoData(allmutations, chrom, start):
    data = {}
    keys = ["chrom","pos","period","prior","family","child","phenotype","posterior","newallele","mutsize","inparents","poocase"]
    for k in keys: data[k] = []
    x = tabix.open(allmutations)
    try:
        records = list(x.query(str(chrom), start-1, start))
    except tabix.TabixError:
        MSG("No mutations found for %s:%s"%(chrom, start))
        return pd.DataFrame(data)
    for r in records:
        if int(r[1]) != start: continue
        for i in range(len(keys)):
            data[keys[i]].append(r[i])
    # Fix types
    df = pd.DataFrame(data)
    df["pos"] = df["pos"].apply(int)
    df["period"] = df["period"].apply(int)
    df["prior"] = df["prior"].apply(float)
    df["phenotype"] = df["phenotype"].apply(int)
    df["posterior"] = df["posterior"].apply(float)
    return df

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--allmutations", help="Indexed file of all mutations", type=str, required=True)
    parser.add_argument("--loci", help="Bed file of all loci to analyze", type=str, required=True)
    parser.add_argument("--annotations", help="Bed file of all annotations to add to output", type=str, required=True)
    parser.add_argument("--out", help="Name of output bed file or stdout", type=str, required=True)
    parser.add_argument("--output-mutations", help="Name of output file for mutations passing filters and thresholds", type=str, required=True)
    parser.add_argument("--filter-both-kids", help="Filter families if a locus was called denovo in both kids", action="store_true")
    parser.add_argument("--max-filtered-families", help="Skip locus if at least this many families are filtered", type=int, default=10000)
    parser.add_argument("--min-children", help="Remove loci if less than this many total unaffected children analyzed", type=int, default=0)
    parser.add_argument("--pthresh", help="Posterior threshold to call something a mutation", type=float, required=True)
    parser.add_argument("--chrom", help="Only look at loci on this chromosome", type=str, default=None)
    parser.add_argument("--profile", help="Profile a couple loci with pygraphviz", action="store_true")
    args = parser.parse_args()

    # Check input
    if not os.path.exists(args.allmutations):
        ERROR(args.allmutations + " does not exist")
    if not os.path.exists(args.allmutations + ".tbi"):
        ERROR("No index found for " + args.allmutations)
    if not os.path.exists(args.loci):
        ERROR(args.loci + " not found")
    if not os.path.exists(args.annotations):
        ERROR(args.annotations + " not found")

    # Load annotations
    ann = pd.read_csv(args.annotations, sep="\t")
    anncols = list(ann.columns[3:])

    # Load loci
    loci = pd.read_csv(args.loci, sep="\t", names=["chrom", "start"])
    if args.chrom is not None:
        loci = loci[loci["chrom"].apply(str) == args.chrom]

    # Get output ready
    if args.out == "stdout":
        outf = sys.stdout
    elif args.profile:
        outf = open(os.devnull, "w")
    else:
        outf = open(args.out, "w")
    if args.output_mutations == "stdout":
        outm = sys.stdout
    elif args.profile:
        outm = open(os.devnull, "w")
    else:
        outm = open(args.output_mutations, "w")

    outcols = ["#chrom", "start", "end", "period", \
               "total_children", "total_mutations", "total_mutation_rate", \
               "affected_children", "affected_mutations", "affected_mutation_rate", \
               "unaffected_children", "unaffected_mutations", "unaffected_mutation_rate", \
               "p-value", "children_with_mutations"] + anncols
    outf.write("\t".join(outcols)+"\n")

    outmcols = ["#chrom","pos", "period", "prior", "family", "child", "phenotype", "posterior", "newallele", "mutsize", "inparents", "poocase"]
    outm.write("\t".join(outmcols)+"\n")
    if args.output_mutations != "stdout":
        outm.close()
        outm = open(args.output_mutations, "a") # so pandas can append to it

    # Process each locus
    for i in range(loci.shape[0]):
        chrom = loci.chrom.values[i]
        start = loci.start.values[i]
        loc_ann = ann[(ann["chrom"]==chrom) & (ann["start"]==start)]
        if loc_ann.shape[0] == 0:
            MSG("No annotation data for %s:%s"%(chrom, start))
            continue
        elif loc_ann.shape[0] > 1:
            MSG("Multiple annotation lines for %s:%s"%(chrom, start))
            continue
        else:
            anndata = []
            for col in anncols: anndata.append(loc_ann[col].values[0])
            end = loc_ann["end"].values[0]
        loc_denovo = ReadDenovoData(args.allmutations, chrom, start)
        denovodata = ProcessLocus(loc_denovo, args.pthresh, outm, outmcols, args.filter_both_kids,
                                  args.max_filtered_families, args.min_children)
        if denovodata is None:
            MSG("Skipping locus %s:%s, no data after filtering or too many families failed filter."%(chrom, start))
            continue
        outdata = [chrom, start, end] + denovodata + anndata
        outf.write("\t".join(map(str, outdata))+"\n")
        outf.flush()

    outf.close()
    outm.close()
if __name__ == "__main__":
    main()
