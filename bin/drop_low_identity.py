#!/usr/bin/env python
from __future__ import division
import argparse
import contextlib
import csv
import logging

import pysam

def parse_blast_contaminants(fp):
    r = csv.DictReader(fp)
    for row in r:
        yield row['query']

def main():
    p = argparse.ArgumentParser()
    p.add_argument('input')
    p.add_argument('output')
    p.add_argument('-i', '--min-identity', default=0.90, help="""Remove reads
            which match reference with < [value] [default: %(default)s]""", type=float)
    p.add_argument('-b', '--blast-contaminants', action='append')
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='[%(name)s] %(message)s')
    logger = logging.getLogger('low_identity')

    blast_contaminants = set()
    for f in a.blast_contaminants:
        with open(f) as fp:
            blast_contaminants |= set(parse_blast_contaminants(fp))

    dropped = 0
    processed = 0
    with contextlib.closing(pysam.Samfile(a.input, 'rb')) as input_bam:
        with contextlib.closing(pysam.Samfile(a.output, 'wb', template=input_bam)) as output_bam:
            for read in input_bam:
                processed += 1

                if read.qname in blast_contaminants:
                    dropped += 1
                    continue

                pct_id = 1.0 - read.opt('NM') / read.alen
                if pct_id < a.min_identity:
                    dropped += 1
                else:
                    output_bam.write(read)
    logger.info('Removed %d/%d [%0.2f%%]', dropped, processed, dropped / processed * 100)
    pysam.index(a.output)
    logger.info("Indexed.")

if __name__ == '__main__':
    main()
