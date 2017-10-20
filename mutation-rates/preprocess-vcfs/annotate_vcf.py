#!/usr/bin/env python
"""
Prepare VCF files with TMRCA field for mutea analysis
"""

import argparse
import os
import sys
import tabix
import tempfile
import vcf

def GetWriter(reader):
    tmpdir = tempfile.mkdtemp(prefix="lobstr.")
    tmpfile = os.path.join(tmpdir, "header.vcf")
    f = open(tmpfile, "w")
    for line in reader._header_lines: f.write(line.strip() + "\n")
    f.write("##FORMAT=<ID=TMRCA,Number=1,Type=Float,Description=\"TMRCA of haplotypes.\">\n")
    f.write("#" + "\t".join(reader._column_headers + reader.samples) + "\n")
    f.close()
    writer = vcf.Writer(sys.stdout, vcf.Reader(open(tmpfile, "rb")))
    return writer

def FuzzyEqual(a, b, buffer=10):
    return abs(a-b) <= buffer

def GetTMRCAData(asdt, chrom, start, end):
    records = asdt.query(chrom, start, end)
    tmrcas = {}
    for record in records:
        if record[0] == chrom and FuzzyEqual(int(record[1]), start) and FuzzyEqual(int(record[2]), end):
            sample = record[4]
            tmrca = float(record[3])
            tmrcas[sample] = tmrca
    return tmrcas

def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--invcf", help="HipSTR or lobSTR VCF input file", required=True, type=str)
    parser.add_argument("--asdt", help="Bed files with chrom, start, end, sample, trmca", required=True, type=str)
    args = parser.parse_args()
    invcf = args.invcf
    asdt = args.asdt

    asdt = tabix.open(asdt)
    reader = vcf.Reader(open(invcf, "rb"))
    writer = GetWriter(reader)

    for record in reader:
        tmrcas = GetTMRCAData(asdt, record.CHROM, record.POS, record.INFO["END"])
        if "TMRCA" not in reader.formats:
            record.add_format("TMRCA")
        samp_fmt = vcf.model.make_calldata_tuple(record.FORMAT.split(':'))
        for fmt in samp_fmt._fields:
            if fmt == "TMRCA":
                entry_type = "Float"
                entry_num = 1
            else:
                entry_type = reader.formats[fmt].type
                entry_num = reader.formats[fmt].num
            samp_fmt._types.append(entry_type)
            samp_fmt._nums.append(entry_num)
        new_samples = []
        for sample in record:
            sampdat = []
            for i in range(len(samp_fmt._fields)):
                key = samp_fmt._fields[i]
                if key != "TMRCA":
                    sampdat.append(sample[key])
                else:
                    tmrca = tmrcas.get(sample.sample, -1)
                    sampdat.append(tmrca)
            call = vcf.model._Call(record, sample.sample, samp_fmt(*sampdat))
            new_samples.append(call)
        record.samples = new_samples
        writer.write_record(record)

if __name__ == "__main__":
    main()
