#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
import sys
from sys import getsizeof
import numpy as np
import argparse

def getKaryotype(fname):
    """ zwraca słownik np: {'chr13': 115169878, ... } """
    data = [i.strip().split() for i in open(fname) if i[:3] == 'chr']
    hs = {}
    for i in data:
        hs[i[6]] = int(i[5])
    return hs

def getSegmentsByChr(data, chromosome, karyotype, resolution, ignoreLong=False, prefix='', interactionLengthCutoff=1.8e7):
    signal = np.zeros(karyotype[chromosome]/resolution, dtype=np.uint16)
    data2 = []
    for i in data:
        if i[0] == chromosome and i[3] == chromosome:
            s = (int(i[1]), int(i[5]), int(i[6]))
            data2.append(s)

    N = len(data2)

    maximum = 0
    warned = False      # sprawdza czy user zostal poinformowany o bliskim przekroczeniu rozmiaru zmiennej
    warned2 = False     # sprawdza czy user zostal poinformowany o wykryciu dlugich inteterakcji
    for k, i in enumerate(data2):
        n = i[2]
        interactionLength = i[1]-i[0]+1
        sys.stdout.write( '\r{:.2f} %, Max: {}'.format(k/N*100, maximum) )
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
                    print "\n[ \033[1;33mWARN\033[1;m ] poziom sygnału przekroczył wartość 95% pojemności typu np.uint16! Dalsze zwiększenie może prowadzić do błędów!"
                    warned = True
    outFname = '{}.{}.signal.np.int16'.format(prefix,chromosome)
    print "\nSaving file {}...".format(outFname)
    f = open(outFname, 'wb')
    signal.tofile(f)

def main():
    parser = argparse.ArgumentParser(description='Zapisuje sygnal chapetowy w postaci pliku binarnego numpy array np.int16')
    parser.add_argument('karyotyp', help='karyotyp w formacie analogicznym do circosa')
    parser.add_argument('dane', help='plik z danymi w formacie bed (klastry ChIA-PET)')
    parser.add_argument('chromosom', help='chromosom. Np chr22')
    parser.add_argument('-p', '--prefix', default='', help='prefix zapisu plik. Może być użyty np do wskazania linii komórkowej. Przykładowa nazwa pliku z prefixem: GM.chr14.signal.np.int16')
    parser.add_argument('-r', '--resolution', type=int, default=10000, help='Im mniej tym lepiej. Rozdzieloczosc. 1 oznacza ze PET-Count jest zliczany dla kazdej pozycji w genomie. Może się liczyć ponad dobę! Rozsądna wartości to: 5-100')
    parser.add_argument('-i', '--ignoreLong', action='store_true', help='Ignore interactions longer than 1.8e7.')
    args = parser.parse_args()

    karyotype = getKaryotype(args.karyotyp)
    data = [ i.strip().split() for i in open(args.dane) ]
    S = getSegmentsByChr(data, args.chromosom, karyotype, args.resolution, args.ignoreLong,  args.prefix ) 
    
   
if __name__ == '__main__':
    main()
