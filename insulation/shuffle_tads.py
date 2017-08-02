#!/usr/bin/env python

import random
import sys

def main(boundaries, chrname):
    # compute deltas from supplied boundaries
    deltas = []
    for l in open(boundaries).readlines():
        a,b = map(int, l.strip().split('\t'))
        deltas.append(b - a)
    # shuffle and write
    random.shuffle(deltas)
    locus = 0
    for d in deltas:
        print('%s\t%d\t%d' % (chrname, locus, locus + d))
        locus += d

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: %s boundaries_file chrname' % sys.argv[0])
    else:
        main(sys.argv[1], sys.argv[2])
