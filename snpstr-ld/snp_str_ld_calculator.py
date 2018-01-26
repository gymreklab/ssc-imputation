#!/usr/bin/env python

"""
Calculate LD between an STR and SNP variant

# Test
./snp_str_ld_calculator.py \
  --str-vcf /storage/s1saini/hipstr_rerun/chr5/hipstr.chr5.with.1kg.filtered.vcf.gz \
  --snp-vcf /storage/resources/datasets/SSC_SNP_v2/shapeit.chr5.with.ref.vcf.gz \
  --str-locus 5:153681115 \
  --snp-locus-rsid rs11740474 \
  --use-info-start

TODO add:
- allelic r2 (return list of ld, rather than single)
- all pairwise snp/str in some window
- --str-vcf2 argument if want to compare imputed vs. ground truth
- remove outliers

"""

import argparse
import sys
import vcf
import scipy.stats

START_BUFFER = 10

def PrintLine(str_locus, snp_locus, ld):
    if ld is None: return
    ld_r, ld_pval = ld
    ld_r2 = ld_r**2
    sys.stdout.write("\t".join(map(str, [str_locus, snp_locus, ld_r2, ld_pval]))+"\n")
    sys.stdout.flush()

def CalcLD(str_reader, snp_reader, str_locus, snp_locus, use_info_start=False, samples=[]):
    sample_to_gts = {} # sample -> {"str", "snp"}
    # Find STR record.
    chrom, start = str_locus.split(":")
    start = int(start)
    str_record = None
    if use_info_start:
        records = str_reader.fetch(chrom, start-START_BUFFER, start+START_BUFFER)
        for r in records:
            if r.INFO["START"] == start:
                str_record = r
                break
        if str_record is None:
            sys.stderr.write("ERROR: couldn't find STR record for %s\n"%start)
            return None
    else:
        records = str_reader.fetch(chrom, start, start+1)
        for r in records:
            str_record = r
            if r is None:
                sys.stderr.write("ERROR: couldn't find STR record for %s\n"%start)
                return None
            break
        if str_record is None or str_record.POS != start:
            sys.stderr.write("ERROR: couldn't find STR record for %s\n"%start)
            return None
    for sample in str_record:
        sample_to_gts[sample.sample] = {"STR": None, "SNP": None}
        if sample["GB"]:
            sample_to_gts[sample.sample]["STR"] = sum(map(int, sample["GB"].split("|")))
    # Find SNP record
    snp_chrom, snp_start = snp_locus.split(":")
    snp_start = int(snp_start)
    records = snp_reader.fetch(snp_chrom, snp_start-1, snp_start)
    snp_record = None
    for r in records:
        snp_record = r
        break
    if snp_record is None:
        sys.stderr.write("Could not find SNP locus %s\n"%snp_locus)
        return None
    if snp_record.POS != snp_start:
        sys.stderr.write("ERROR: couldn't find SNP record for %s\n"%snp_start)
        return None
    for sample in snp_record:
        if sample.sample not in sample_to_gts.keys():
            continue
        if sample["GT"]:
            sample_to_gts[sample.sample]["SNP"] = sum(map(int, sample.gt_alleles))
    str_data = []
    snp_data = []
    for sample in sample_to_gts:
        if len(samples)>0 and sample not in samples: continue
        if sample_to_gts[sample]["STR"] is not None and sample_to_gts[sample]["SNP"] is not None:
            str_data.append(sample_to_gts[sample]["STR"])
            snp_data.append(sample_to_gts[sample]["SNP"])
    return scipy.stats.pearsonr(str_data, snp_data)

def FindSnpLocus(snp_reader, str_locus, snp_locus_rsid, snp_region_start, snp_region_end):
    chrom, start = str_locus.split(":")
    start = int(start)
    try:
        records = snp_reader.fetch(chrom, snp_region_start, snp_region_end)
    except:
        sys.stderr.write("Could not fetch records for %s\n"%snp_locus_rsid)
        return None
    for r in records:
        if r.ID == snp_locus_rsid or str(r.POS) in snp_locus_rsid:
            return "%s:%s"%(r.CHROM, r.POS)
    return None

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--str-vcf", help="VCF with STR genotypes", type=str, required=True)
    parser.add_argument("--snp-vcf", help="VCF with SNP genotype", type=str, required=True)
    parser.add_argument("--str-locus", help="chr:start of STR locus", type=str, required=False)
    parser.add_argument("--snp-locus", help="chr:start of SNP locus", type=str, required=False)
    parser.add_argument("--snp-locus-rsid", help="rsid of SNP locus", type=str, required=False)
    parser.add_argument("--loci-file", help="File with chr,start(STR),start(SNP)", type=str, required=False)
    parser.add_argument("--loci-file-rsid", help="File with chr,start(STR),rsid. Optional columns snp_start, snp_end", type=str, required=False)
    parser.add_argument("--pairwise-snpstr", help="Calculate pairwise LD for all SNPs within maxdist of the STR", action="store_true")
    parser.add_argument("--max-dist", help="Don't consider snp/str more than this many bp apart", type=int, default=100000)
    parser.add_argument("--use-info-start", help="Match STR start on INFO/START (not POS)", action="store_true")
    parser.add_argument("--samples", help="Only consider samples in this file", type=str, required=False)
    args = parser.parse_args()

    # Open readers
    str_reader = vcf.Reader(open(args.str_vcf, "rb"))
    snp_reader = vcf.Reader(open(args.snp_vcf, "rb"))

    # Get samples
    samples = []
    if args.samples is not None:
        samples = [item.strip() for item in open(args.samples, "r").readlines()]

    if args.str_locus:
        str_locus = args.str_locus
        snp_locus = None
        if args.snp_locus:
            snp_locus = args.snp_locus
            snpid = snp_locus
            ld = CalcLD(str_reader, snp_reader, str_locus, snp_locus, \
                        use_info_start=args.use_info_start, samples=samples)
        elif args.snp_locus_rsid:
            snpid = args.snp_locus_rsid
            str_start = int(str_locus.split(":")[1])
            snp_region_start = str_start - args.max_dist
            snp_region_end = str_start + args.max_dist
            snp_locus = FindSnpLocus(snp_reader, args.str_locus, args.snp_locus_rsid, \
                                     snp_region_start, snp_region_end)
            ld = CalcLD(str_reader, snp_reader, args.str_locus, snp_locus, \
                        use_info_start=args.use_info_start, samples=samples)
            if snp_locus is None:
                sys.stderr.write("ERROR: Couldn't find SNP locus within %s of %s\n"%(args.max_dist, args.str_locus))
        else:
            sys.stderr.write("ERROR: No SNP locus specified. Use --snp-locus or --snp-locus-rsid\n")
        PrintLine(str_locus, snpid, ld)

    # Keep track of snp id positions in case we see same one many times
    snp_rsid_to_pos = {}
    if args.loci_file_rsid is not None or args.loci_file is not None:
        has_rsid = False
        fname = args.loci_file
        if args.loci_file_rsid is not None:
            fname = args.loci_file_rsid
            has_rsid = True
        with open(fname, "r") as f:
            for line in f:
                items = line.strip().split()
                chrom = items[0]
                str_start = int(items[1])
                str_locus = "%s:%s"%(chrom, str_start)
                snp_loci = []
                if has_rsid:
                    rsids = items[2].split(",")
                    try:
                        snp_region_start = int(items[3])
                    except IndexError:
                        snp_region_start = str_start-max_dist
                    try:
                        snp_region_end = int(items[4])
                    except IndexError:
                        snp_region_end = str_start+max_dist
                    for rsid in rsids:
                        snp_locus = None
                        try:
                            snp_locus = snp_rsid_to_pos[rsid]
                            if snp_locus is None: continue
                        except KeyError:
                            snp_locus = FindSnpLocus(snp_reader, str_locus, rsid, \
                                                     snp_region_start, snp_region_end)
                            snp_rsid_to_pos[rsid] = snp_locus
                        if snp_locus is None:
                            sys.stderr.write("Could not find %s\n"%rsid)
                            snp_rsid_to_pos[rsid] = None
                            continue
                        snp_loci.append(snp_locus)
                else:
                    snp_start = int(items[2])
                    snp_loci = ["%s:%s"%(chrom, item) for item in snp_start.split(",")]
                for i in range(len(snp_loci)):
                    snp_locus = snp_loci[i]
                    snpid = snp_locus
                    if args.loci_file_rsid is not None:
                        snpid = rsids[i]
                    ld = CalcLD(str_reader, snp_reader, str_locus, snp_locus, \
                                use_info_start=args.use_info_start, \
                                samples=samples)
                    PrintLine(str_locus, snpid, ld)

    if args.pairwise_snpstr:
        str_reader2 = vcf.Reader(open(args.str_vcf, "rb"))
        for str_record in str_reader2:
            str_locus = "%s:%s"%(str_record.CHROM, str_record.POS)
            region_start = max([0, str_record.POS - args.max_dist])
            region_end = str_record.POS + args.max_dist
            try:
                snp_records = snp_reader.fetch(str_locus.split(":")[0], region_start, region_end)
            except ValueError:
                sys.stderr.write("ERROR fetching SNP records for STR locus %s\n"%str_locus)
                continue
            snp_loci = ["%s:%s"%(record.CHROM, record.POS) for record in snp_records]
            for snp_locus in snp_loci:
                ld = CalcLD(str_reader, snp_reader, str_locus, snp_locus, \
                            samples=samples)
                PrintLine(str_locus, snp_locus, ld)
main()
