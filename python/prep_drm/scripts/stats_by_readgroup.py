#!/usr/bin/env python
"""
Classify mutations for read groups in a BAM file.
"""

import argparse
import collections
import contextlib
import csv
import logging
import sys

import numpy as np
import pysam
from scipy.stats.mstats import mquantiles


def main():
    p = argparse.ArgumentParser()
    p.add_argument('bamfile', type=pysam.Samfile)
    p.add_argument('-o', '--output', type=argparse.FileType('w'),
                   default=sys.stdout)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    alens, nms = collections.defaultdict(list), collections.defaultdict(list)
    with contextlib.closing(a.bamfile), a.output as fp:
        for read in a.bamfile:
            rg = read.opt('RG')
            nms[rg].append(read.opt('NM'))
            alens[rg].append(read.alen)

        sam_header = a.bamfile.header
        rg_infos = {i['ID']: i for i in sam_header['RG']}

        w = csv.writer(fp, lineterminator='\n')
        w.writerow(('sample', 'library', 'n', 'mean_aligned_length',
            'mean_mismatches', 'mean_pct_id', 'median_pct_id',
            'pct_id_2.5', 'pct_id_97.5'))
        for rg in sorted(alens):
            rg_info = rg_infos[rg]
            alen = np.array(alens[rg])
            nm = np.array(nms[rg])
            pct_id = (alen - nm).astype(np.float64) / alen
            row = [rg, rg_info['LB'], len(alen), np.mean(alen), np.mean(nm), np.mean(pct_id), np.median(pct_id)]
            row.extend(mquantiles(pct_id, [0.025, 0.975]))
            w.writerow(row)

if __name__ == '__main__':
    main()
