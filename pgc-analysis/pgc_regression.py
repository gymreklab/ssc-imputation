#!/usr/bin/env python
"""
Regression analysis between STRs and SCZ phenotype in PGC
"""

import argparse
import vcf

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--strvcf", help="VCF file with STR genotypes", required=True, type=str)
    parser.add_argument("--sampledata", help="Tab file with pheno/covars for each sample", required=True, type=str)
    parser.add_argument("--out", help="Prefix for output files. Default: stdout", required=False, type=str)
    parser.add_argument("--region", help="Restrict analysis to STRs in this region.", type=str)
    args = parser.parse_args()

    str_reader = 
if __name__ == "__main__":
    main()
