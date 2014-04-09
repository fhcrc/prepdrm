#!/usr/bin/env python
from __future__ import division

import argparse
import collections
from contextlib import closing
import functools
import logging
import math
import sys

import pysam

PHRED_OFFSET = 33

def keep_position(q, qlen, is_reverse):
    if is_reverse:
        q = qlen - q
    return q <= 275

def main():
    p = argparse.ArgumentParser()
    p.add_argument('bam_file', metavar='bam!', type=pysam.Samfile)
    p.add_argument('proxy_reference', help="""Make consensus for samples
                   aligning to this reference""")
    p.add_argument('-m', '--min-coverage', default=20, type=int)
    p.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'))
    p.add_argument('-n', '--name')
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO)

    counterdict = functools.partial(collections.defaultdict, collections.Counter)
    result = collections.defaultdict(counterdict)

    with closing(a.bam_file) as bam:
        h = bam.header
        sample_map = {i['ID']: i['SM'] for i in h['RG']}
        is_reverse = {i['ID']: 'd2' in i['ID'] for i in h['RG']}
        reference_map = {i: n for i, n in enumerate(bam.references)}

        if a.proxy_reference not in bam.references:
            raise KeyError('reference {0} not found in {1}'.format(
                           a.proxy_reference, bam.references))


        for read in bam:
            ref_name = reference_map[read.tid]
            if ref_name != a.proxy_reference:
                continue
            rg = read.opt('RG')
            is_rev = is_reverse[rg]
            sample = sample_map[rg]
            qbases = read.query
            qlen = read.qlen
            for q, r in read.aligned_pairs:
                if q and r and keep_position(q, qlen, is_rev):
                    result[sample][r][qbases[q]] += 1
    for sample, counts in result.iteritems():
        min_ref = min(counts)
        max_ref = max(counts)
        logging.info('%s: %d-%d', sample, min_ref, max_ref)
        bases = []
        qualities = []

        for i in xrange(min_ref, max_ref + 1):
            c = counts[i]
            if not c:
                bases.append('N')
                qualities.append(chr(1 + PHRED_OFFSET))
                continue
            most_prevalent = max(c, key=c.get)
            bases.append(most_prevalent)
            p_wrong = 1.0 - c[most_prevalent] / float(sum(c.values()))
            if p_wrong == 0.0:
                qualities.append(chr(40 + PHRED_OFFSET))
            else:
                qualities.append(chr(int(-10 * math.log(p_wrong, 10)) + PHRED_OFFSET))

        a.output.write('@{0}\n{1}\n+{0}\n{2}\n'.format(sample + '_cons', ''.join(bases), ''.join(qualities)))

if __name__ == '__main__':
    main()
