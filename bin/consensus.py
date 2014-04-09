#!/usr/bin/env python
"""
Quick consensus given two BAM files
"""
import argparse
import contextlib
import collections
import sys

import pysam

FORWARD = (919, 1186)
REVERSE = (1186, 1442)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('forward', type=pysam.Samfile)
    p.add_argument('reverse', type=pysam.Samfile)
    p.add_argument('-r', '--reference', default='MB2059_pol')
    p.add_argument('-n', '--name', default='Consensus')
    a = p.parse_args()

    seq = []
    for (start, end), bam in [(FORWARD, a.forward), (REVERSE, a.reverse)]:
        with contextlib.closing(bam):
            for pileup in bam.pileup(a.reference, start, end, max_depth=1e6):
                if pileup.pos < start or pileup.pos >= end:
                    continue
                c = collections.Counter()
                for p in pileup.pileups:
                    if p.indel:
                        continue
                    c[p.alignment.seq[p.qpos]] += 1
                seq.append(max(c, key=c.get))
    print '>{0}\n{1}'.format(a.name, ''.join(seq))


if __name__ == '__main__':
    main()
