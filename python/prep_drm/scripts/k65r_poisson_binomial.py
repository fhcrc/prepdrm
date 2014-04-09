#!/usr/bin/env python
import argparse
import collections

import numpy as np
import pysam

POS = 954

def p_all_failed(probs):
    return sum(np.log(probs))

def decode_phred_quality(char, offset=33):
    return ord(char) - offset

def main():
    p = argparse.ArgumentParser()
    p.add_argument('bamfile', type=pysam.Samfile)
    a = p.parse_args()

    piled = a.bamfile.pileup('MB2059_pol', start=POS, max_depth=1e6)

    pos, count, pileups = [(i.pos, i.n, i.pileups)
                           for i in piled
                           if i.pos == POS][0]

    c = collections.Counter()
    qualities = collections.defaultdict(list)

    for pu_read in pileups:
        alignment = pu_read.alignment
        base = alignment.query[pu_read.qpos]
        phred_quality = alignment.qqual[pu_read.qpos]
        q = decode_phred_quality(phred_quality)
        c[base] += 1
        qualities[base].append(q)

    g_quals = np.array(qualities['G'], np.float64)
    g_probs = 10. ** (- g_quals / 10.0)
    log_p_failed = p_all_failed(g_probs)
    print c
    print "(log,raw) probability all failed: ", log_p_failed, np.exp(log_p_failed)

if __name__ == '__main__':
    main()
