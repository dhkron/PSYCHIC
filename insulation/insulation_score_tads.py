#!/usr/bin/env python

import sys
import os

H_FORMAT = 'bin%d|tt|%s:%d-%d'

def gen_new_matrix(matrix, chrname, res):
    n = len(matrix)
    header = '\t'.join([H_FORMAT % (i, chrname, 1 + i * res, 1 + (i + 1) * res) \
                        for i in range(n)])
    new_matrix = '\t' + header + '\n'
    for i, row in enumerate(matrix):
        new_matrix += H_FORMAT % (i, chrname, 1 + i * res, 1 + (i + 1) * res)
        new_matrix += '\t' + row
    return new_matrix

def main(matrix_file, chrname, res):
    # add header to input matrix 
    matrix = open(matrix_file).readlines()
    fname = matrix_file + '.header' 
    boundaries_fname = os.path.basename(matrix_file) + '.is520001.ids320001.insulation.boundaries'
    o = gen_new_matrix(matrix, chrname, res)
    f = open(fname, 'w')
    f.write(o)
    f.close()

    # run TAD calling, and convert result to our format
    if not os.path.exists(boundaries_fname):
        os.system('./matrix2insulation.pl -i %s' % fname)
    boundaries = open(boundaries_fname).readlines()[1:]
    boundaries_bins = [int(x.split('\t')[4]) for x in boundaries]


    # clean
    #if raw_input('clean? ') == 'y': 
    os.system('rm -f %s.is520001.ids320001.insulation*' % os.path.basename(matrix_file))

    for i in range(1, len(boundaries_bins)):
        print('%s\t%d\t%d' % (chrname, boundaries_bins[i-1] * res, boundaries_bins[i] * res))



if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.stderr.write("Usage: %s input_matrix.txt chrname resolution" % sys.argv[0])
    else:
        main(sys.argv[1], sys.argv[2], int(sys.argv[3]))


