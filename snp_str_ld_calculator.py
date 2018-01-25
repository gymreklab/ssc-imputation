#!/usr/bin/env python

"""
Calculate LD between an STR and SNP variant

# Test
./snp_str_ld_calculator.py \
  --str-vcf /storage/s1saini/hipstr_rerun/chr5/hipstr.chr5.with.1kg.filtered.vcf.gz \
  --snp-vcf /storage/resources/datasets/SSC_SNP_v2/shapeit.chr5.with.ref.vcf.gz \
  --str-locus 5:153681115 \
  --snp-locus rs11740474
"""

import argparse
import sys
import vcf

def CalcLD(str_reader, snp_reader, str_locus, snp_locus):
    return 0 # TODO

def FindSnpLocus(str_locus, snp_locus_rsid, max_dist):
    return 0 # TODO

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--str-vcf", help="VCF with STR genotypes", type=str, required=True)
    parser.add_argument("--snp-vcf", help="VCF with SNP genotype", type=str, required=True)
    parser.add_argument("--str-locus", help="chr:start of STR locus", type=str, required=False)
    parser.add_argument("--snp-locus", help="chr:start of SNP locus", type=str, required=False)
    parser.add_argument("--snp-locus-rsid", help="rsid of SNP locus", type=str, required=False)
    parser.add_argument("--loci-file", help="File with chr,start(STR),start(SNP)", type=str, required=False)
    parser.add_argument("--loci-file-rsid", help="File with chr,start(STR),rsid", type=str, required=False)
    parser.add_argument("--max-dist", help="Don't consider snp/str more than this many bp apart", type=int, default=1000000)
    args = parser.parse_args()

    # Open readers
    str_reader = vcf.Reader(open(args.str_vcf, "rb"))
    snp_reader = vcf.Reader(open(args.snp_vcf, "rb"))

    if args.str_locus:
        str_locus = args.str_locus
        snp_locus = None
        if args.snp_locus:
            snp_locus = args.snp_locus
            ld = CalcLD(str_reader, snp_reader, str_locus, snp_locus)
        elif args.snp_locus_rsid:
            snp_locus = FindSnpLocus(args.str_locus, args.snp_locus_rsid, args.max_dist)
            ld = CalcLD(str_reader, snp_reader, args.str_locus, snp_locus)
            if snp_loc is None:
                sys.stderr.write("ERROR: Couldn't find SNP locus within %s of %s\n"%(args.max_dist, args.str_locus))
        else:
            sys.stderr.write("ERROR: No SNP locus specified. Use --snp-locus or --snp-locus-rsid\n")

        sys.stdout.write("\t".join(map(str, [str_locus, snp_locus, ld]))+"\n")

    # TODO read from file

main()
