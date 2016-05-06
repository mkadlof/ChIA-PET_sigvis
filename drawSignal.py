#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import sys
from matplotlib import pyplot as plt
import numpy as np
import argparse

def getKaryotype(fname):
    """returns dictionary e.g.: {'chr13': 115169878, ... } """
    data = [i.strip().split() for i in open(fname) if i[:3] == 'chr']
    hs = {}
    for i in data:
        hs[i[6]] = int(i[5])
    return hs

def drawSignals(karyotype, chromosome, signalFiles, segments):
    for f in signalFiles:
        y = np.fromfile(f, dtype=np.uint16)
        if y.size > 1e6:
            every = int(y.size//1e6)
            print "[ \033[1;33mWARN\033[1;m ] Your data from {} are very big so they were probed every {:d} point so they could fit in circa 1 000 000 points. Sorry :(".format(f.name, every)
            y = y[0::every]
        x = np.linspace(0,karyotype[chromosome],len(y))
        print "Plotting {}".format(f.name)
        bla = plt.plot(x,y, '-', label=f.name)
        color = bla[-1].get_color()
        plt.fill_between(x,y, color = color, alpha=0.1)

    if segments:
        tmp = [i.strip().split() for i in segments ]
        segByPrzemek = [int(i[1]) for i in tmp if i[0] == chromosome]
        plt.plot(segByPrzemek,[0 for i in xrange(len(segByPrzemek))],'d',color='yellow',linestyle="none" )

    plt.legend()
    plt.grid()
    plt.title("ChIA-PET signal for {}".format(chromosome))
    plt.show()
    
def main():
    parser = argparse.ArgumentParser(description='Script for plotting ChIA-PET signal')
    parser.add_argument('karyotype', help='karyotype in Circos format')
    parser.add_argument('chromosome', help='chromosome. Np chr22')
    parser.add_argument('signalFiles', type=argparse.FileType('r'), nargs='+', help="filenames with signal saved in numpy np.uint16 binary format")
    parser.add_argument('-s', '--segments', type=argparse.FileType('r'), help="name of file with segments")
    args = parser.parse_args()
    
    karyotype = getKaryotype(args.karyotype)
    drawSignals(karyotype, args.chromosome, args.signalFiles, args.segments)
    
if __name__ == '__main__':
    main()
