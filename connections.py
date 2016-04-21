#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import sys
from sys import getsizeof
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import patches
from math import log, floor, sqrt
from scipy import ndimage
import argparse

def getKaryotype(fname):
    """returns dictionary e.g.: {'chr13': 115169878, ... } """
    data = [i.strip().split() for i in open(fname) if i[:3] == 'chr']
    hs = {}
    for i in data:
        hs[i[6]] = int(i[5])
    return hs

def drawLinks(data, karyotype, chromsome, ax):
    for i in data:
        archWidth = i[5] - i[1]
        arcCenter_x = archWidth/2+i[1]
        arcHeight = i[6]                # PET-Count
        arc = patches.Arc((arcCenter_x,0), archWidth, arcHeight)
        ax.add_patch(arc)


def main():
    parser = argparse.ArgumentParser(description="Plots connections diagram")
    parser.add_argument('-s', '--segments', type=argparse.FileType('r'), help="filename with segmentami")
    parser.add_argument('karyotype', help='karyotype in Circos format')
    parser.add_argument("clusters", help="file with ChIA-PET clusters")
    parser.add_argument("chromosome", help="np chr22" )
    args = parser.parse_args()

    chromosome = args.chromosome
    print "Openning karyotype..."
    karyotype = getKaryotype(args.karyotype)
    print "Read clusters..."
    clusters = [ i.strip().split() for i in open(args.clusters) ]
    clusters = [ [i[0], int(i[1]), int(i[2]), i[3], int(i[4]), int(i[5]), int(i[6]) ] for i in clusters]
    clusters = [ i for i in clusters if i[0] == chromosome and i[3] == chromosome]

    maxPetCount = max(zip(*clusters)[6])
    
    print "Plotting..."
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.set_xlim(0, karyotype[chromosome])    
    ax1.set_ylim(0, maxPetCount*0.55)
    drawLinks(clusters, karyotype, chromosome, ax1)

    segments = args.segments
    if segments:
        tmp = [i.strip().split() for i in segments ]
        segments = [int(i[1]) for i in tmp if i[0] == chromosome]
        ax1.plot(segments,[0 for i in xrange(len(segments))],'d',color='yellow',linestyle="none" )
        for i in segments:
            ax1.axvline(x=i, alpha=0.3, color="green")

    plt.show()
    

if __name__ == '__main__':
    main()
