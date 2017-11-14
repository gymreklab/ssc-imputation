#!/usr/bin/env python

SFARIFILE="/storage/mgymrek/ssc-denovos/denovos2/other-data/SFARI-Gene_cnvs_export13-11-2017.csv"
CYTOFILE="/storage/mgymrek/ssc-denovos/denovos2/other-data/cytoBand.txt.gz"
OUTFILE="/storage/mgymrek/ssc-denovos/denovos2/other-data/SFARI-CNV-annotated.tab"

import pandas as pd
import sys

sfari = pd.read_csv(SFARIFILE)
cyto = pd.read_csv(CYTOFILE, sep="\t", names=["chrom","start","end","band","np"])
cyto["name"] = cyto.apply(lambda x: x["chrom"][3:]+x["band"], 1)

def GetStartEnd(x, cyto):
    loc = (x["cnv-locus"])
    locX = loc.replace("p","?")
    locX = locX.replace("q","?")
    chrom = locX.split("?")[0]
    if "-" in loc:
        locstart = loc.split("-")[0]
        locend = chrom + loc.split("-")[1]
    else: 
        locstart = loc
        locend = loc
    try:
        start = cyto[cyto["name"].apply(lambda x: locstart in x)]["start"].values[0]
    except:
        sys.stderr.write("Could not convert %s\n"%locstart)
        return "NA"
    try:
        end = cyto[cyto["name"].apply(lambda x: locend in x)]["end"].values[-1]
    except:
        sys.stderr.write("Could not convert %s\n"%locend)
        return "NA"
    return "%s:%s-%s"%(chrom, start, end)

sfari["loc"] = sfari.apply(lambda x: GetStartEnd(x, cyto), 1)
sfari = sfari[sfari["loc"] != "NA"]

sfari["chrom"] = sfari["loc"].apply(lambda x: x.split(":")[0])
sfari["start"] = sfari["loc"].apply(lambda x: x.split("-")[0].split(":")[1])
sfari["end"] = sfari["loc"].apply(lambda x: x.split("-")[1])

sfari[["chrom","start","end","number-case-individuals"]].to_csv(OUTFILE, sep="\t", index=False)
