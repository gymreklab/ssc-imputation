#!/usr/bin/env python

"""
Scale mutation rates to have same mean as calibration set (e.g. CODIS)
"""

import argparse
import numpy as np
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--infile", help="Mutea estimates", type=str, required=True)
    parser.add_argument("--outfile", help="Output adjusted estimates here", type=str, required=True)
    parser.add_argument("--scale", help="Scale for mutation rates. Output will be MUTEA/scale", type=float, required=True)
    parser.add_argument("--gamma", help="Scale factor for standard errors", type=float, required=True)
    args = parser.parse_args()

    # Load data
    columns = ["chrom","start","end","est_logmu_ml","est_beta_ml","est_pgeom_ml", \
                   "est_logmu_stderr","numsamples_ml","center"]

    data = pd.read_csv(args.infile, sep="\t", names=columns)
    sys.stderr.write("before filter %s\n"%data.shape[0])

    # Set stderr nan to -1 so we know who they are
    data.ix[np.isnan(data["est_logmu_stderr"]), "est_logmu_stderr"] = -1 # set to -1, will filter later
    data.ix[np.isinf(data["est_logmu_stderr"]), "est_logmu_stderr"] = -1 # set to -1, will filter later
    sys.stderr.write("after filter stderr %s\n"%data.shape[0])

    # Scale mutation rates using CODIS scaling factor
    data["est_logmu_ml"] = data["est_logmu_ml"] - np.log10(args.scale)

    # Scale stderrs
    data["est_logmu_stderr"] = data.apply(lambda x: x.est_logmu_stderr*abs(x.est_logmu_ml)*args.gamma, 1)

    # Output results
    data[columns].to_csv(args.outfile, sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()
