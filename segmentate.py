#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import numpy as np
import argparse
from os.path import basename, dirname

def getKaryotype(fname):
    """returns dictionary e.g.: {'chr13': 115169878, ... } """
    data = [i.strip().split() for i in open(fname) if i[:3] == 'chr']
    hs = {}
    for i in data:
        hs[i[6]] = int(i[5])
    return hs

def performSegmentation(signal, karyotype, chromosome, CUTOFF=3e5, CUTOFF2=9e5):
    """signal - numpy.uint16 binary file
       karyotype - description of chromosomes, same as Circos file format
       chromosome - symbol from last column of karyotype file. e.g. 'chr22'
       CUTOFF - minimal segment size
       CUTOFF2 - closest segment shouldn't be closer than this value
    """

    directory = dirname(signal)
    infname = basename(signal)
    
    signal = open(signal, 'r')

    karyotype = getKaryotype(karyotype)
    signal = np.fromfile(signal, dtype=np.uint16)
    if signal.size > 1e6:
        every = int(signal.size/1e6)
        print "[ \033[1;33mWARN\033[1;m ] Your data from {} are very big so they were probed every {:d} point so they could fit in circa 1 000 000 points. Sorry :(".format(f.name, every)
        signal = signal[0::every]
    
    segmentsCoords = []

    inSegment = False
    for i in xrange(len(signal)):
        if inSegment and signal[i] == 0:
            continue
        elif not inSegment and signal[i] == 0:
            inSegment = True
            begin = i
        elif inSegment and signal[i] != 0:
            inSegment = False
            end = i
            segmentsCoords.append([begin, end])
        elif not inSegment and signal[i]:
            continue
    if signal[-1] == 0:
        end = len(signal)
        segmentsCoords.append([begin, end])

    # usuniecie poczatku i konca
    del segmentsCoords[0]
    del segmentsCoords[-1]

    # przeskalowanie
    n = signal.size
    dlugoscChromosomu = karyotype[chromosome]
    skala = dlugoscChromosomu/n
    for i in xrange(len(segmentsCoords)):
        for j in xrange(2):
            segmentsCoords[i][j] *= skala

    segments = []
    for i in segmentsCoords:
        begin = i[0]
        end = i[1]
        length = end - begin
        position = (length/2)+begin
        if length > CUTOFF:
            segments.append([(begin,end), position, length])
   
    segments2 = []
    for i in xrange(len(segments)-1):
        if segments[i][1] not in segments2:
            if segments[i+1][1]-segments[i][1] > CUTOFF2: # jeśli dwa segmenty są zbyt blisko siebie
                if segments[i+1][2] > segments[i][2]:
                    segments2.append(segments[i+1][1])
                else:
                    segments2.append(segments[i][1])
            else:
                segments2.append(segments[i][1])

    outFname = "{}.segments.txt".format(chromosome)

    if directory != '': outPath = "{}/{}".format(directory, outFname)
    else: outPath = "{}".format(outFname)

    f = open(outPath,'w')
    for i in segments2:
        f.write('{} {} {}\n'.format(chromosome, i, i ))
    print "{} segments identified.".format(len(segments2))
    print '{} saved.'.format(outPath)
    
def main():
    parser = argparse.ArgumentParser(description='Perform ChIA-PET signal segmentation')
    parser.add_argument('karyotype', help='karyotype in Circos format')
    parser.add_argument('signalFile', help="ChIA-PET clusters signal in  np.uint16 binary format")
    parser.add_argument('chromosome', help='chromosme. e.g. chr22')
    parser.add_argument('-c1', '--cutoff1', type=int, default=3e5,  help='minimal segment size. Default: 3e5 bp')
    parser.add_argument('-c2', '--cutoff2', type=int, default=9e5,  help='closest segment shouldn\'t be closer than this value. Default: 9e5 bp')
    args = parser.parse_args()

    performSegmentation(args.signalFile, args.karyotype, args.chromosome, CUTOFF=args.cutoff1, CUTOFF2=args.cutoff2)

if __name__ == '__main__':
    main()
