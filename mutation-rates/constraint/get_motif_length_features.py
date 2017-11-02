#!/usr/bin/env python
"""
get motif and length features for each item in hipstr reference
"""

import sys
import pyfaidx
import operator

try:
    hipreffile = sys.argv[1]
    reffafile = sys.argv[2]
except:
    sys.stderr.write(__doc__)
    sys.exit(1)

ERROR = "N"
nucToNumber={"A":0,"C":1,"G":2,"T":3}

def CheckMotif(motif):
    # If all consist of the same nucleotide, error
    if len(motif) > 1 and len(''.join(set(motif))) == 1: return False
    return True

def GetRepseq(repseq):
    """ Get canonical STR sequence, considering both strands """
    repseq_f = getCanonicalMS(repseq)
    repseq_r = getCanonicalMS(reverseComplement(repseq))
    repseq = compareString(repseq_f, repseq_r)
    return repseq

def getCanonicalMS(repseq):
    """ Get canonical STR sequence """
    size = len(repseq)
    canonical = repseq
    for i in range(size):
        newseq = repseq[size-i:]+repseq[0:size-i]
        for j in range(size):
            if nucToNumber[newseq[j]] < nucToNumber[canonical[j]]:
                canonical = newseq
            elif nucToNumber[newseq[j]] > nucToNumber[canonical[j]]:
                break
    return canonical

def compareString(seq1,seq2):
    """ Compare two strings alphabetically """
    size = len(seq1)
    for i in range(size):
        if nucToNumber[seq1[i]] < nucToNumber[seq2[i]]:
            return seq1
        if nucToNumber[seq1[i]] > nucToNumber[seq2[i]]:
            return seq2
    return seq1

def reverseComplement(seq):
    """ Get the reverse complement of a nucleotide string """
    newseq = ""
    size = len(seq)
    for i in range(len(seq)):
        char = seq[len(seq)-i-1]
        if char == "A": newseq += "T"
        if char == "G": newseq += "C"
        if char == "C": newseq += "G"
        if char == "T": newseq += "A"
    return newseq

def GetMostCommonSeq(seq, period):
    seq = seq.upper()
    if len(seq) < period: return ERROR
    motifs = {}
    for i in range(len(seq)-period):
        subseq = seq[i:i+period]
        if CheckMotif(subseq):
            motifs[subseq] = motifs.get(subseq, 0) + 1
    if len(motifs.keys()) == 0: return ERROR
    bestmotif = max(motifs.iteritems(), key=operator.itemgetter(1))[0]
    return GetRepseq(bestmotif)

def GetUninterruptedLength(chrom, start, end, motif, genome):
    seq = genome[chrom][start-1:end].upper()
    rmotif = reverseComplement(motif)
    sub1 = GetUninterruptedLength_(seq, motif)
    sub2 = GetUninterruptedLength_(seq, rmotif)
    return max([sub1, sub2])

def GetUninterruptedLength_(seq, motif):
    loc = 0
    starts = []
    ends = []
    ind = seq.find(motif, loc) # get first one
    if ind != -1:
        starts.append(ind)
        loc = ind + len(motif)
    while True:
        ind = seq.find(motif, loc)
        if ind == -1: # didn't find anymore, end run
            if len(starts) > 0:
                ends.append(loc-1)
            break
        if ind > loc: # imperfect, end run and start new one
            ends.append(min([loc-1, len(seq)-1]))
            starts.append(ind)
        loc = ind + len(motif)
    intervals = zip(starts, ends)
    intervals = map(list, intervals)
    # For each interval, try to extend
    maxlen = 0
    bestint = [-1, -1]
    for i in range(len(intervals)):
        subseq = seq[intervals[i][0]:intervals[i][1]+1]
        prefix = seq[max([intervals[i][0]-len(motif), 0]):intervals[i][0]]
        suffix = seq[intervals[i][1]+1:min([intervals[i][1]+1+len(motif), len(seq)])]
        length = intervals[i][1]-intervals[i][0] # full match
        for j in range(1, len(motif)+1): # match prefix
            if len(prefix) < j: break
            if prefix[-1*j] == motif[-1*j]:
                intervals[i][0] -= 1
            else:
                break
        for j in range(len(motif)): # match suffix
            if len(suffix) <= j: break
            if suffix[j] == motif[j]:
                intervals[i][1] += 1
            else:
                break
        if intervals[i][1]-intervals[i][0]+1 > maxlen:
            maxlen = intervals[i][1]-intervals[i][0]+1
            bestint = intervals[i]
    return maxlen

# Load reference genome
genome = pyfaidx.Fasta(reffafile, as_raw=True)

with open(hipreffile, "r") as f:
    for line in f:
        items = line.strip().split() 
        chrom = items[0]
        start = int(items[1])
        end = int(items[2])
        period = int(items[3])
        seq = str(genome[chrom][start:end+1])
        motif = GetMostCommonSeq(seq, period)
        length = end-start+1
        if "N" not in motif:
            uninterrupted_length = GetUninterruptedLength(chrom, start, end, motif, genome)
        else:
            uninterrupted_length = 0
        sys.stdout.write("\t".join(map(str, [chrom, start, end, motif, length, uninterrupted_length]))+"\n")
