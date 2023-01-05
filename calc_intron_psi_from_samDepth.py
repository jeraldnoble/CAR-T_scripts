#!/usr/bin/env python3
import sys,re,io
from collections import defaultdict
import numpy as np

def main():
    depth, sj, bed = getArgs()
    depth = storeDepth(depth = depth)
    sj, sample = storeSJ(sj = sj)
    bed = storeBed(bed = bed)
    calcPsi(depth = depth, sj = sj, bed = bed, sample = sample)

def getArgs():
    import argparse
    parser = argparse.ArgumentParser(description = "calc psi of introns")
    parser.add_argument('-depth', action = 'store', type = str, help = 'samtools depth output')	
    parser.add_argument('-sj', action = 'store', type = str, help = 'star sj out file')	
    parser.add_argument('-bed', action = 'store', type = str, help = 'bed file of intron coords')	
    args = parser.parse_args()
    depth, sj, bed = args.depth, args.sj, args.bed
    return(depth, sj, bed)

def storeBed(bed):
    bedDict = defaultdict(dict)
    for line in open(bed, 'r'):
        arr = line.split('\t')
        bedDict[arr[0]][(arr[1] + '-' + arr[2])] = arr[3]
    return(bedDict)

def storeDepth(depth):
    depthDict = defaultdict(dict)
    for line in open(depth, 'r'):
        line = line.strip()
        arr = line.split('\t')
        depthDict[arr[0]][arr[1]] =int(arr[2])
    return(depthDict)

def storeSJ(sj):
    sample = sj.split('_filt')[0]
    sjDict = defaultdict(dict)
    for line in open(sj, 'r'):
        line = line.strip()
        arr = line.split('\t')
        sjDict[arr[0]][str(arr[1] + '-' + arr[2])] = int(arr[3])
    return(sjDict, sample)
    

def calcPsi(depth, bed, sj, sample):
    ## get sj reads 28936201 28936512
    ## exon boundaries 28936153 28936201 28936512 28936600
    events = []
    psiVals = []
    for chromo in bed:
        for coords in bed[chromo]:
            event = chromo + ':' + coords
            events.append(event)
            strand = bed[chromo][coords]
            left, right = coords.split('-')
            left, right = int(left), int(right)
            left = left + 1
            if (str(left) + '-' + str(right)) in sj[chromo]:
                sjCov = sj[chromo][str(left) + '-' + str(right)]
                irCov = []
                for i in range(left, right+1):
                    irCov.append(depth[chromo][str(i)])
                irCov = np.median(irCov)
                psi = irCov / (irCov + sjCov)
                psi = round(psi, 3)
                psiVals.append(str(psi))
            else:
                psiVals.append("NA")
    print("sample\t" + '\t'.join(events))
    print(sample + '\t' + '\t'.join(psiVals))

if __name__ == "__main__":
	main()
