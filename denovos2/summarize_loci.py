#!/usr/bin/env python
"""
Using filtered denovo calls, output stats per locus
"""

import argparse
import pandas as pd
import tabix

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--allmutations", help="Indexed file of all mutations", type=str, required=True)
    parser.add_argument("--loci", help="Bed file of all loci to analyze", type=str, required=True)
    parser.add_argument("--annotations", help="Bed file of all annotations to add to output", type=str, required=True)
    parser.add_argument("--out", help="Name of output bed file", type=str, required=True)
    args = parser.parse_args()

if __name__ == "__main__":
    main()
