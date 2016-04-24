#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import sys
from sys import getsizeof
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

def getSegmentsByChr(data, chromosmee, karyotype, resolution=10000, ignoreLong=False, prefix='', interactionLengthCutoff=4e5):

    directory = dirname(data)
    infname = basename(data)

    data = [ i.strip().split() for i in open(data) ]
    signal = np.zeros(karyotype[chromosmee]/resolution, dtype=np.uint16)
    data2 = []
    for i in data:
        if i[0] == chromosmee and i[3] == chromosmee:
            s = (int(i[1]), int(i[5]), int(i[6]))
            data2.append(s)

    N = len(data2)

    maximum = 0
    warned = False      # sprawdza czy user zostal poinformowany o bliskim przekroczeniu rozmiaru zmiennej
    warned2 = False     # sprawdza czy user zostal poinformowany o wykryciu dlugich inteterakcji
    for k, i in enumerate(data2):
        n = i[2]
        interactionLength = i[1]-i[0]+1
        sys.stdout.write( '\r{:.2f} %, Max signal value: {}'.format(k/N*100, maximum) )
        sys.stdout.flush()

        # Filtrowanie po długości interackcji
        if interactionLength > interactionLengthCutoff:
            if ignoreLong:
                if not warned2:
                    print "\n[ \033[1;33mWARN\033[1;m ] Very long interaction detected! They will be ignored."
                    warned2 = True
                continue
            else:
                if not warned2: 
                    print "\n[ \033[1;33mWARN\033[1;m ] Very long interaction detected! They will be included."
                    warned2 = True

        for j in xrange(i[0], i[1], resolution):
            j = int(j/resolution)
            signal[j] += n
            if signal[j] > maximum:
                maximum = signal[j]
                if maximum > 62258:
                    print "\n[ \033[1;33mWARN\033[1;m ] level of maximal signal value exceded 95% of maximal np.uint16 capacity! There is a risk of overflow!"
                    warned = True
    
    if prefix != '': outFname = '{}.{}.signal.np.uint16'.format(prefix,chromosmee)
    else: outFname = '{}.signal.np.uint16'.format(chromosmee)
    
    if directory != '': outPath = "{}/{}".format(directory, outFname)
    else: outPath = "{}".format(outFname)
    
    print "\nSaving file {}...".format(outFname)
    f = open(outPath, 'wb')
    signal.tofile(f)

def main():
    parser = argparse.ArgumentParser(description='Saved ChIA-PET signal into numpy uint16 one-dimension binary array format in any desired resolution.')
    parser.add_argument('karyotype', help='karyotype in Circos format')
    parser.add_argument('data', help='file with data (ChIA-PET clusters)')
    parser.add_argument('chromosome', help='chromosme. e.g. chr22')
    parser.add_argument('-p', '--prefix', default='', help='prefix of output file name. e.g.: GM.chr14.signal.np.uint16 (\'GM\' is prefix here)' )
    parser.add_argument('-r', '--resolution', type=int, default=10000, help='The lower, the better. Resolution equal 1 mean that PET-Count is counted for each position in genome. It may last over 24h! Resonable values ar between 5 and 10000')
    parser.add_argument('-i', '--ignoreLong', action='store_true', help='Ignore interactions longer than 18 Mbases.')
    args = parser.parse_args()

    karyotype = getKaryotype(args.karyotype)
    print "Processing {}...".format(args.chromosome)
    S = getSegmentsByChr(args.data, args.chromosome, karyotype, args.resolution, args.ignoreLong,  args.prefix ) 
    
   
if __name__ == '__main__':
    main()
