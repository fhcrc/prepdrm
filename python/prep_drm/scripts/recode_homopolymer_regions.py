"""
Recode homopolymer regions, adding/subtracting bases if the difference is
within some threshold
"""
import argparse
from contextlib import closing
import itertools
import logging
import operator

import pysam

from .. import gff3, pwalign, rle, samutil
from ..util import isolate_region

log = logging.getLogger(__name__.rpartition('.')[-1])
BAM_CMATCH = 0

class MissingReferenceSequenceError(ValueError):
    pass

def parse_gff3(fp):
    """
    Read mutations from a GFF3 file
    """
    records = gff3.parse(fp)
    records = (record for record in records
               if record.type == 'possible_base_call_error')
    return records

def rle_bases(rle_encoded):
    return ''.join(i.c for i in rle_encoded)

def correct_region(aligned_pairs, max_diff=1):
    aligned_pairs = list(aligned_pairs)

    # Check for fully-aligning sequence
    if all(i.qpos is not None and i.rpos is not None for i in aligned_pairs):
        return aligned_pairs

    ref_bases = ''.join(i.rbase for i in aligned_pairs if i.rbase)
    read_bases = ''.join(i.qbase for i in aligned_pairs if i.qbase)

    if not read_bases:
        return aligned_pairs

    ref_rle = rle.encode(ref_bases)
    ref_rle_bases = rle_bases(ref_rle)
    read_rle = rle.encode(read_bases)
    read_rle_bases = rle_bases(read_rle)

    aref, aqry, _ = pwalign.pw_global(ref_rle_bases, read_rle_bases)

    if '-' in aref or '-' in aqry:
        # No perfect alignment after RLE
        return aligned_pairs

    # Maximum number of hp changes to make
    max_hp = max(abs(r.length - q.length) for q, r in zip(read_rle, ref_rle))
    if max_hp == 0 or max_hp > max_diff:
        # Too much to correct / nothing to correct
        return aligned_pairs

    # Okay, actually fix
    ref_idx = itertools.count(next(i.rpos for i in aligned_pairs if i.rpos is not None)).next
    qry_idx = itertools.count(next(i.qpos for i in aligned_pairs if i.qpos is not None)).next

    qqual = {i.qpos: i.qual for i in aligned_pairs if i.qpos is not None}

    result = [samutil.AlignedPair(qpos=qry_idx(),
                                  rpos=ref_idx(),
                                  qbase=q.c,
                                  rbase=r.c,
                                  qual=chr(33+5),
                                  cigar_op=BAM_CMATCH)
              for q, r in zip(read_rle, ref_rle)
              for i in xrange(r.length)]

    result = [i._replace(qual=qqual.get(i.qpos, chr(33+1)))
              if i.qpos is not None else i
              for i in result]
    return result


def hp_correct(read, reference, regions, max_diff=1):
    all_pairs = list(samutil.all_pairs_iter(read, reference))
    ap = [(i.qpos, i.rpos) for i in all_pairs]
    in_region = frozenset(i for (s, e) in regions
                          for i in isolate_region(ap, s, e))
    grouped = itertools.groupby(all_pairs, lambda x: (x.qpos, x.rpos) in in_region)
    result = [i for g, v in grouped
              for i in (correct_region(v, max_diff=max_diff) if g else v)]

    read.seq = ''.join(i.qbase for i in result if i.qbase is not None)
    read.qual = ''.join(i.qual for i in result if i.qbase is not None)
    read.cigar = [(op, sum(True for i in v))
                  for op, v in itertools.groupby(i.cigar_op for i in result)]
    return read


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-m', '--max-diff', type=int, default=1, help="""Maximum
                   difference between observed HP length and expected HP length
                   to correct""")
    p.add_argument('reference', type=pysam.Fastafile)
    p.add_argument('gff3_file', type=argparse.FileType('r'))
    p.add_argument('input_bam', type=pysam.Samfile)
    p.add_argument('output_bam')
    a = p.parse_args()

    with a.gff3_file as fp:
        k = operator.attrgetter('seqid')
        regions = {s: list(g)
                   for s, g in itertools.groupby(sorted(parse_gff3(fp), key=k), k)}

    with closing(a.input_bam), closing(a.reference), \
            closing(pysam.Samfile(a.output_bam, 'wb', template=a.input_bam)) as out_bam:
        available_references = {i: (r, a.reference.fetch(r)) for i, r in enumerate(a.input_bam.references)}

        for read in a.input_bam:
            ref_name, ref_bases = available_references[read.tid]
            if not ref_bases:
                raise MissingReferenceSequenceError(ref_name)
            reg = [(region.start0, region.end) for region in regions[ref_name]]
            read = hp_correct(read, ref_bases, reg, max_diff=a.max_diff)
            out_bam.write(read)
