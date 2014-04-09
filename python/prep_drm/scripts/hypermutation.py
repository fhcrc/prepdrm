"""
Looking for evidence of additional G->A mutations in reads classified as M184I
"""
import argparse
import collections
import contextlib
import csv
import logging
import sys

import pysam

from .classify_mutations import parse_gff3, classify_read
from ..util import opener, window

log = logging.getLogger(__name__.rpartition('.')[-1])


def list_gg_to_ag(read, reference, mutation, min_quality=20):
    query = read.query
    qual = [ord(char) - 33 for char in read.qqual]

    counts = collections.Counter()
    dimers = list(window(read.aligned_pairs, 2))
    rcounts = collections.Counter()
    for (q1, r1), (q2, r2) in dimers:
        if any(i is None for i in (q1, q2, r1, r2)):
            continue
        overlaps_mut = any(r >= mutation.codon_start0 and
                           r < mutation.codon_start0 + 3
                           for r in (r1, r2))
        low_qual = qual[q1] < min_quality or qual[q2] < min_quality
        if overlaps_mut or low_qual:
            continue

        r = reference[r1], reference[r2]
        q = query[q1], query[q2]
        rcounts[r] += 1

        if r == ('G', 'G'):
            counts[q] += 1

    total = sum(counts.itervalues())

    a = counts.get(('A', 'G'), 0)
    return {'ref_gg': total,
            'qry_gg': counts.get(('G', 'G'), 0),
            'qry_ag': a,
            'gg_prop_mutated': float(a) / total if total else None}


def list_g_to_a(read, reference, mutation, min_quality=20):
    query = read.query
    qual = [ord(char) - 33 for char in read.qqual]

    counts = collections.Counter(query[q] for q, r in read.aligned_pairs
                                 if q is not None and
                                 r is not None and
                                 reference[r] == 'G' and
                                 qual[q] >= min_quality and
                                 (r < mutation.codon_start0 or r > mutation.codon_start0 + 2))

    total = sum(counts.itervalues())
    a = counts.get('A', 0)
    return {'ref_g': total,
            'qry_g': counts.get('G', 0),
            'qry_a': a,
            'g_prop_mutated': float(a) / total if total else None}

def main():
    p = argparse.ArgumentParser()
    p.add_argument('fasta', type=pysam.Fastafile)
    p.add_argument('gff3_file', type=argparse.FileType('r'))
    p.add_argument('bamfile', type=pysam.Samfile)
    p.add_argument('-q', '--min-quality', default=20, type=int, help="""Minimum
                   quality score to consider base [default: %(default)d]""")
    p.add_argument('-o', '--output', type=opener('w'),
                   default=sys.stdout)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s: %(message)s')

    with a.gff3_file as fp:
        mutations = {i.reference: i for i in parse_gff3(fp) if i.name == 'M184VI'}

    with contextlib.closing(a.bamfile), a.output as fp:
        headers = ['read', 'read_group', 'amino_acid', 'classification',
                   'ref_g', 'qry_g', 'qry_a', 'g_prop_mutated',
                   'ref_gg', 'qry_gg', 'qry_ag', 'gg_prop_mutated']
        w = csv.DictWriter(fp, headers, lineterminator='\n')
        w.writeheader()

        references = [a.fasta.fetch(reference_name) for reference_name in a.bamfile.references]
        mutations = [mutations[reference_name] for reference_name in a.bamfile.references]

        for read in a.bamfile:
            reference_name = a.bamfile.references[read.tid]
            reference_sequence = references[read.tid]
            m = mutations[read.tid]
            c = classify_read(read, m)
            if c is None:
                # No classification
                continue
            hm = list_g_to_a(read, reference_sequence, m, a.min_quality)
            result = {'read': read.qname}
            result.update({k: v for k, v in c._asdict().iteritems()
                           if k in {'read_group', 'amino_acid', 'classification'}})
            result.update(hm)
            result.update(list_gg_to_ag(read, reference_sequence, m, a.min_quality))

            w.writerow(result)


if __name__ == '__main__':
    main()
