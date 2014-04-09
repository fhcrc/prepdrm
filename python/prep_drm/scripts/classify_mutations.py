#!/usr/bin/env python
"""
Classify mutations for read groups in a BAM file,
using mutations identified by a GFF3 file
"""

import argparse
import collections
import contextlib
import csv
import itertools
import re
import logging
import sys

from Bio.Data import CodonTable
import pysam

from .. import gff3
from ..util import isolate_region

log = logging.getLogger(__name__.rpartition('.')[-1])

CODON_TABLE = CodonTable.ambiguous_dna_by_name['Standard'].forward_table

_MUT_FIELDS = ('parent', 'mut')
_ML_FIELDS = ('reference', 'name', 'codon_start1', 'wt', 'mutations')

MUTATION_PATTERN = re.compile(r'([A-Z]|\+\d+)')


class UnrecognizedMutationError(ValueError):
    pass


class MutationLocation(collections.namedtuple('MutationLocation', _ML_FIELDS)):
    def classify_amino_acid(self, amino_acid, read, codon_indels):
        # If any of the mutations at this location are positive,
        # it's a mutant.
        if any(m.is_mutant(amino_acid, read)
               for m in self.mutations
               if not codon_indels or m.supports_indels):
            return 'mut'
        # Indels / invalid AA -> non-coding
        elif amino_acid in '-!X' or codon_indels:
            return 'non_coding'
        elif amino_acid in self.wt:
            return 'wt'
        return 'other'

    @property
    def codon_start0(self):
        """0-based codon start position"""
        return self.codon_start1 - 1


class Mutation(collections.namedtuple('Mutation', _MUT_FIELDS)):
    """
    Description of a mutation
    """
    is_insertion = False
    supports_indels = False

    def is_mutant(self, amino_acid, read):
        if amino_acid in self.mut:
            return True

class Insertion(Mutation):
    """
    Insertion relative to reference
    """
    is_insertion = True
    supports_indels = True

    def __init__(self, *args, **kwargs):
        super(Insertion, self).__init__(*args, **kwargs)
        if not self.mut.startswith('+'):
            raise UnrecognizedMutationError(self.mut)

    @property
    def length(self):
        return int(self.mut[1:])

    def is_mutant(self, amino_acid, read):
        loc = self.parent.codon_start0 + 2

        ap = read.aligned_pairs
        it = itertools.dropwhile(lambda (q, r): r < loc, iter(ap))
        try:
            q, r = next(it)
        except StopIteration:
            return
        if not (q and r):
            return

        insertions = list(itertools.takewhile(lambda (q, r): r is None, it))
        return len(insertions) == 3 * self.length


def parse_gff3(fp):
    """
    Read mutations from a GFF3 file
    """
    records = gff3.parse(fp)

    p = re.compile(r'^([A-Z]+)(?:\d+)((?:[A-Z]|\+\d)+)$')
    for record in records:
        # Multiple nucleotide polymorphisms only (SO:0001013)
        if record.type != 'MNP':
            continue
        assert record.end == record.start + 2
        attributes = record.attribute_dict()
        name = attributes['Name']
        m = p.match(name)
        if not m:
            raise UnrecognizedMutationError("Unrecognized name: " + name)
        wt, mut = m.groups()

        mloc = MutationLocation(record.seqid, name, record.start, wt, [])
        for m in MUTATION_PATTERN.findall(mut):
            muts = ((Insertion if m.startswith('+') else Mutation)(mloc, m)
                    for m in MUTATION_PATTERN.findall(mut))
        mloc.mutations.extend(muts)
        yield mloc

ClassificationResult = collections.namedtuple('ClassificationResult',
                                              ('read_group', 'amino_acid',
                                               'classification', 'sequence',
                                               'codon_sequence', 'any_indels',
                                               'whole_codon_indels'))


def classify_read(read, mutation, context=3, allow_indels_in_context=False):
    """
    Classify a read.

    :param AlignedRead read: AlignedRead, from pysam
    :param Mutation mutation: Mutation description
    :param bool context: Number of nucleotides on each side to consider
    :param bool allow_indels_in_context: Should indels in the context regions
    cause a read to be classified as non-coding?
    """
    codon_start = mutation.codon_start0

    try:
        rg = read.opt('RG')
    except KeyError:
        rg = None

    # Find position covered
    region_start = codon_start - context
    region_end = codon_start + 2 + context + 1
    region_aligned_pairs = isolate_region(read.aligned_pairs,
                                          region_start,
                                          region_end)
    if not region_aligned_pairs:
        return


    n_leading_del = region_aligned_pairs[0][1] - region_start
    assert n_leading_del >= 0
    if n_leading_del:
        region_aligned_pairs = [(None, region_start + i)
                                for i in xrange(n_leading_del)] + region_aligned_pairs

    s = []
    codon_str = []
    query = read.query
    for q, r in region_aligned_pairs:
        q_nt = query[q] if q is not None else '-'
        if r is None:
            q_nt = q_nt.lower()
        s.append(q_nt)

        codon_nt = q_nt
        if r is not None and (r < codon_start or r > codon_start + 2):
            codon_nt = codon_nt.lower()
        codon_str.append(codon_nt)

    codon = ''.join(query[q] if q is not None else '-'
                    for q, _ in isolate_region(region_aligned_pairs,
                                               codon_start, codon_start + 3))

    if all(i == '-' for i in s) or len(codon) < 3:
        return None

    any_indels = any(i is None or j is None for i, j in region_aligned_pairs)

    g = lambda (i, j): (i is None, j is None)
    segment_lengths = [sum(True for i in v)
                       for _, v in itertools.groupby(region_aligned_pairs, g)]
    non_codon_indels = any(l % 3 != 0 for l in segment_lengths)

    # translate
    amino_acid = '!'
    if (not non_codon_indels or allow_indels_in_context) and len(codon) % 3 == 0 and len(codon) >= 3:
        try:
            amino_acid = CODON_TABLE.get(codon[:3], 'X')
        except CodonTable.TranslationError:
            amino_acid = 'X'

    classification = mutation.classify_amino_acid(amino_acid, read, any_indels and not non_codon_indels)

    return ClassificationResult(rg, amino_acid, classification, ''.join(s),
                                ''.join(codon_str), any_indels,
                                any_indels and not non_codon_indels)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('bamfile', type=pysam.Samfile)
    p.add_argument('gff3_file', type=argparse.FileType('r'))
    p.add_argument('-o', '--output', type=argparse.FileType('w'),
                   default=sys.stdout)
    p.add_argument('--context', type=int, default=3, help="""Number of
            nucleotides on each side of mutation location to show""")
    p.add_argument('--allow-indels-in-context', default=False,
            action='store_true', help="""Should indels in the context region
            lead to a non_coding classification?""")
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s:%(name)s: %(message)s')

    with a.gff3_file as fp:
        mutations = list(parse_gff3(fp))

    available_references = set(a.bamfile.references)

    with contextlib.closing(a.bamfile), a.output as fp:
        w = csv.writer(fp, lineterminator='\n')
        w.writerow(('location', 'sample', 'amino_acid', 'classification',
                    'sequence', 'codon_highlighted', 'any_indels',
                    'whole_codon_indels', 'count'))
        # Fetch
        for i, mutation in enumerate(mutations):
            if mutation.reference not in available_references:
                log.warn("Skipping %s: %s not in BAM references", mutation.name, mutation.reference)
                continue
            log.info('%s: %s [%d/%d]', mutation.reference, mutation.name, i+1, len(mutations))
            reads = a.bamfile.fetch(mutation.reference, mutation.codon_start0,
                                    mutation.codon_start1)

            counts = collections.Counter()
            for read in reads:
                r = classify_read(read, mutation, context=a.context,
                                  allow_indels_in_context=a.allow_indels_in_context)
                if r:
                    counts[r] += 1

            log.info('%s: coverage at %s (%s) = %s',
                     mutation.reference,
                     mutation.name,
                     mutation.codon_start1,
                     sum(counts.values()))
            items = sorted(counts.iteritems(),
                           key=lambda (i, c): (i.read_group, -c))
            rows = ((mutation.name, i.read_group, i.amino_acid,
                     i.classification, i.sequence, i.codon_sequence,
                     'yes' if i.any_indels else 'no',
                     'yes' if i.whole_codon_indels else 'no', count)
                    for i, count in items)
            w.writerows(rows)


if __name__ == '__main__':
    main()
