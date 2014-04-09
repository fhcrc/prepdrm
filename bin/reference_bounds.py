#!/usr/bin/env python
"""
Generates a JSON file with start and end locations of the PrEP 454 amplicon
translated to each sequence.  If a sequence does not cover the entire amplicon,
the maximum covered position is used instead.
"""

import argparse
import operator
import logging
import json
import sys

from Bio import SeqIO

from prep_drm.scripts.translate_gff3 import convert_coords

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--begin', type=int, default=914)
    p.add_argument('--end', type=int, default=1514)
    p.add_argument('-r', '--reference', default='MB2059_pol')
    p.add_argument('fasta', type=argparse.FileType('r'))
    p.add_argument('-o', '--outfile', type=argparse.FileType('w'), default=sys.stdout)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    with a.fasta as fp:
        sequences = {i.id: str(i.seq) for i in SeqIO.parse(fp, 'fasta')}

    ref = sequences[a.reference]
    result = {}
    second = operator.itemgetter(1)
    for k, v in sequences.iteritems():
        if k == a.reference:
            result[k] = (a.begin, a.end)
        else:
            # ref, query coord
            coords = list(convert_coords(ref, v, False))
            begin = next((j for i, j in coords if i == a.begin), None)
            if begin is None:
                ref_begin, begin = min(coords, key=second)
                begin = begin or 0
                logging.warn("%s: Using minimum value %s -> %s", k, ref_begin, begin)

            end = next((j for i, j in coords if i == a.end), None)
            if end is None:
                ref_end, end = max(coords, key=second)
                logging.warn("%s: Using maximum value %s -> %s", k, ref_end, end)

            result[k] = (begin, end)

    with a.outfile as ofp:
        json.dump(result, ofp, indent=2)
        ofp.write('\n')


if __name__ == '__main__':
    main()
