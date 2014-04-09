#!/usr/bin/env python
"""
Given an alignment containing a reference sequence, and the name of a second
sequence, translate the coordinates of features to the second sequence
"""
import argparse
import csv
import logging
import sys
import itertools

from Bio import AlignIO

from .. import gff3


# Field indices
GFF_REF_IDX = 0
GFF_START_IDX = 3
GFF_END_IDX = 4
GFF_INFO_IDX = 8

GFF_CSV_RULES = {'quoting': csv.QUOTE_NONE,
                 'lineterminator': '\n',
                 'delimiter': '\t'}

def convert_coords(ref, query, one_based=True):
    """
    Convert coordinates between a reference and query sequence

    By default, convert to 1-based indices (GFF3 style)

    :param ref: String reference
    :param query: String query
    """
    base = 1 if one_based else 0
    if len(ref) != len(query):
        raise ValueError('Lengths do not match')
    ref_idx = itertools.count(base)
    query_idx = itertools.count(base)

    for ref_char, query_char in itertools.izip_longest(ref, query):
        yield (next(ref_idx) if ref_char != '-' else None,
               next(query_idx) if query_char != '-' else None)

class UnableToTranslateError(ValueError):
    pass

def translate_gff3_feature(feature_row, query, sequence_map):
    """
    Translate a single GFF3 feature to ``query``'s coordinates
    """
    gff_row = feature_row
    reference = sequence_map[gff_row.seqid]
    translation = dict(convert_coords(str(reference.seq), str(query.seq)))
    start_idx, end_idx = gff_row.start, gff_row.end

    t_start = translation[start_idx]
    t_end = translation[end_idx]

    if not t_start:
        raise UnableToTranslateError(
                'Reference start {0} has no equivalent'.format(start_idx))
    if not t_end:
        raise UnableToTranslateError(
                'Reference end {0} has no equivalent'.format(end_idx))

    attr = gff_row.attribute_dict()
    gff_row = gff_row._replace(seqid=query.id, start=t_start, end=t_end)\
            .update_attributes(ID='_'.join((attr['ID'], query.id)))

    logging.info("[%s] Translated (%d,%d) to (%d,%d)\n%s - %s\n%s - %s",
            gff_row.attributes, start_idx, end_idx, t_start, t_end,
            str(reference.seq).translate(None, '-')[start_idx-1:end_idx], reference.id,
            str(query.seq).translate(None, '-')[t_start-1:t_end], query.id)

    return gff_row

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('alignment', type=argparse.FileType('r'))
    p.add_argument('gff3', type=argparse.FileType('r'))
    p.add_argument('-o', '--outfile', type=argparse.FileType('w'),
            default=sys.stdout, help="""Path to output [default: stdout]""")
    p.add_argument('-f', '--alignment-format', default='fasta',
            help="""Format of input alignment [default: %(default)s]""")
    p.add_argument('-q', '--query', help="""Name of query sequence [default: all]""")
    p.add_argument('-k', '--keep-original', help="""Keep original records [default: %(default)s""",
            action='store_true', default=False)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    with a.alignment:
        alignment = AlignIO.read(a.alignment, a.alignment_format)

    sequences = {sequence.id: sequence for sequence in alignment}

    with a.outfile as ofp, a.gff3 as gff_fp:
        w = csv.writer(ofp, **GFF_CSV_RULES)

        # Write comment lines unchaged
        for g, v in itertools.groupby(gff_fp, lambda line: line.startswith('#')):
            if g:
                # Comment
                for line in v:
                    ofp.write(line)
            else:
                # Translate non-comment lines
                rows = gff3.parse(v)

                for row in rows:
                    if a.keep_original:
                        w.writerow(row)
                    if a.query:
                        non_query = [sequences[a.query]]
                    else:
                        non_query = [s for s in sequences.values() if s.id != row[0]]
                    for query in non_query:
                        try:
                            w.writerow(translate_gff3_feature(row, query, sequences))
                        except UnableToTranslateError as e:
                            logging.warn(e)

if __name__ == '__main__':
    main()
