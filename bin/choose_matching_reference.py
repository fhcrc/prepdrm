#!/usr/bin/env python
import argparse
import csv
import itertools
import re
import sys

from prep_drm import gff3

import pysam

def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))

def main():
    p = argparse.ArgumentParser()
    p.add_argument('fasta', type=pysam.Fastafile)
    p.add_argument('gff3', type=argparse.FileType('r'))
    p.add_argument('-p', '--pattern', default=r'.*MB2059.*',
                   help="""Pattern to search""")
    p.add_argument('-o', '--output', default=sys.stdout, type=argparse.FileType('w'))

    a = p.parse_args()

    with a.gff3:
        gff_records = [r._asdict() for r in gff3.parse(a.gff3)
                       if r.type == 'possible_base_call_error' and
                          r.attribute_dict()['Name'] == 'K65R Homopolymer context']

    pat = re.compile(a.pattern)
    def not_candidate(x):
        return pat.match(x['seqid']) is None
    records = sorted(gff_records, key=not_candidate)

    for record in records:
        record['sequence'] = a.fasta.fetch(reference=record['seqid'], start=record['start']-1, end=record['end'])
        assert record['sequence'], record
    a.fasta.close()

    candidates, queries = [list(v) for _, v in itertools.groupby(records, not_candidate)]

    header = ['sequence', 'min_dist'] + [i['seqid'] for i in candidates]

    with a.output as fp:
        w = csv.writer(fp, lineterminator='\n')
        w.writerow(header)
        for q in records:
            distances = [hamming(q['sequence'], r['sequence']) for r in candidates]
            min_i = min(xrange(len(distances)), key=lambda i: distances[i])
            row = [q['seqid'], candidates[min_i]['seqid']] + distances
            w.writerow(row)


if __name__ == '__main__':
    main()
