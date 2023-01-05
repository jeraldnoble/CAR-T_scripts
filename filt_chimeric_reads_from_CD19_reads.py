#!/usr/bin/env python3
import sys,re,io
from collections import defaultdict

def main():
    sam, bed = getArgs()
    bed = storeBed(bed = bed)
    filtSam(sam = sam, bed = bed)

def getArgs():
    import argparse
    parser = argparse.ArgumentParser(description = "remove potential chimeric reads from bam file")
    parser.add_argument('-sam', action = 'store', type = str, help = 'sam file for CD19 loci')	
    parser.add_argument('-bed', action = 'store', type = str, help = 'bed file of parent gene loci for intron retention events')	
    args = parser.parse_args()
    sam, bed = args.sam, args.bed
    return(sam, bed)

def storeBed(bed):
    bedDict = defaultdict(dict)
    fh = open(bed, 'r')
    for line in fh:
        line = line.strip()
        arr = line.split()
        if arr[0] in bedDict:
            bedDict[arr[0]].append((arr[1] + ':' + arr[2]))
        else:
            bedDict[arr[0]] = (arr[1] + ':' + arr[2])
    return(bedDict)

def filtSam(sam, bed):
    for line in open(sam, 'r'):
        line = line.strip()
        if line.startswith("@"):
            print(line)
            continue
        arr = line.split('\t')
        left = int(bed[arr[2]].split(':')[0])
        right = int(bed[arr[2]].split(':')[1])
        span = (right - left) + 2000
        start = int(arr[3])
        cigar = re.split(r'[a-zA-Z]',arr[5])
        readLen = 0
        for i in cigar:
            if i != '':
                readLen += int(i)
        readEnd = left + readLen
        if int(arr[3]) > (left - 2000) and readEnd < (right + 2000):
                print(line)

if __name__ == "__main__":
	main()
