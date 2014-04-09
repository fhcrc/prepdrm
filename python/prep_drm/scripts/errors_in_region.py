#!/usr/bin/env python
"""
Check errors by position
"""

import argparse
import collections
import contextlib
import csv
import logging
import sys

import pysam

from .. import gff3
from ..util import isolate_region

POSS_BASE_CALL_ERROR = 'possible_base_call_error' # SO:0000701

def parse_gff3(fp):
    """
    Read mutations from a GFF3 file
    """
    records = gff3.parse(fp)
    return [i for i in records if i.type == POSS_BASE_CALL_ERROR]

def count_errors(read, feature, reference_base_map):
    try:
        rg = read.opt('RG')
    except KeyError:
        rg = None

    region_aligned_pairs = isolate_region(read.aligned_pairs, feature.start - 1, feature.end)
    substitutions = 0
    del_count = sum(True for q, _ in region_aligned_pairs if q is None)
    ins_count = sum(True for _, r in region_aligned_pairs if r is None)
    for q, r in region_aligned_pairs:
        if q is None or r is None:
            continue
        q_nt = read.query[q] if q is not None else '-'
        r_nt = reference_base_map.get(r, '-')
        if q_nt != r_nt:
            substitutions += 1

    return rg, substitutions, ins_count, del_count


def main():
    p = argparse.ArgumentParser()
    p.add_argument('ref_fasta', type=pysam.Fastafile)
    p.add_argument('bamfile', type=pysam.Samfile)
    p.add_argument('gff3_file', type=argparse.FileType('r'))
    p.add_argument('-o', '--output', type=argparse.FileType('w'),
                   default=sys.stdout)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    with a.gff3_file as fp:
        features = parse_gff3(fp)

    available_references = set(a.bamfile.references)

    with contextlib.closing(a.bamfile), \
         a.output as fp, \
         contextlib.closing(a.ref_fasta):
        w = csv.writer(fp, lineterminator='\n')
        w.writerow(('sample', 'reference', 'feature', 'ref_sequence',
                    'position', 'ref_base', 'n_aligned', 'n_del', 'n_match',
                    'n_mismatch', 'n_ins'))
        # Fetch
        for feature in features:
            feature_id = feature.attribute_dict()['ID']
            feature_name = feature.attribute_dict()['Name']
            if feature.seqid not in available_references:
                logging.warn('Skipping %s: %s not available', feature.seqid, feature_id)
                continue
            start = feature.start - 1
            end = feature.end
            reference_bases = a.ref_fasta.fetch(feature.seqid, start, end)
            reference_base_map = {start + i: b for i, b in enumerate(reference_bases)}

            logging.info('%s %s - %s', feature.seqid, feature_name, reference_bases)

            pileups = a.bamfile.pileup(feature.seqid, start=feature.start - 1, end=feature.end - 1,
                                       max_depth=1e6, fastafile=a.ref_fasta)

            for pileup in pileups:
                if pileup.pos < feature.start - 1 or pileup.pos >= feature.end:
                    continue
                ref_base = reference_base_map[pileup.pos]

                r = collections.defaultdict(collections.Counter)
                for read in pileup.pileups:
                    d = r[read.alignment.opt('RG')]
                    d['n_aligned'] += 1
                    if read.indel > 0:
                        d['n_ins'] += 1
                    if read.is_del:
                        d['n_del'] += 1
                        continue
                    base = read.alignment.seq[read.qpos]
                    if base == ref_base:
                        d['n_match'] += 1
                    else:
                        d['n_mismatch'] += 1
                for sample, d in r.iteritems():
                    w.writerow([sample, feature.seqid, feature_name, reference_bases,
                                pileup.pos - feature.start + 1, ref_base,
                                d['n_aligned'], d['n_del'], d['n_match'], d['n_mismatch'], d['n_ins']])

if __name__ == '__main__':
    main()
