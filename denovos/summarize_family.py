#!/usr/bin/env python
"""
Summarize family de novo results

TODO: filters?
"""

import argparse
import pandas as pd
import sys

def ProcessLine(x, outf, prob, family):
    probs = [x["post_nomut"],x["post_mut_affected"],x["post_mut_unaffected"]]
    if max(probs) < prob:
        sys.stderr.write("Skiping %s %s %s\n"%(family, x["chrom"], x["pos"]))
        return
    items = map(str, [family, x["chrom"], x["pos"]])
    if probs[0] == max(probs):
        outf.write("\t".join(items + ["0", "0", "0", "0"])+"\n")
    elif probs[1] == max(probs):
        if x["mend_lengths"]:
            outf.write("\t".join(items + ["1", "0", "0", "0"])+"\n")
        else:
            outf.write("\t".join(items + ["0", "0", "1", "0"])+"\n")
    else:
        if x["mend_lengths"]:
            outf.write("\t".join(items + ["0", "1", "0", "0"])+"\n")
        else:
            outf.write("\t".join(items + ["0", "0", "0", "1"])+"\n")
    return


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--infile", help="Input file from process_denovo.sh", type=str, required=True)
    parser.add_argument("--outfile", help="Output summary file name", type=str, required=True)
    parser.add_argument("--prob", help="Posterior prob threshold", type=float, required=True)
    parser.add_argument("--family", help="Family name", type=str, required=False, default="NA")
    args = parser.parse_args()

    data = pd.read_csv(args.infile, sep="\t")
    
    outcols = ["family", "chrom", "pos", "mut_affected_length", "mut_unaffected_length", \
               "mut_affected_seq","mut_unaffected_seq"]
    f = open(args.outfile, "w")
    f.write("\t".join(outcols)+"\n")
    data.apply(lambda x: ProcessLine(x, f, args.prob, args.family), 1)
    f.close()

if __name__ == "__main__":
    main()
