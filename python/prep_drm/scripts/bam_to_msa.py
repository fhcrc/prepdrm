"""
Convert a BAM file to a multiple sequence alignment, removing insertions
"""
import argparse
import contextlib
import logging
import sys

import pysam

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--begin', type=int, help="""Start index""")
    p.add_argument('--end', type=int, help="""End index""")
    p.add_argument('--min-reads', type=int, help="""If --begin and --end are
            not specified, find begin and and using first and last position
            with at least MIN_READS aligned bases. [default: %(default)d]""", default=1)
    p.add_argument('-r', '--reference', required=True)
    p.add_argument('reference_fasta', type=pysam.Fastafile)
    p.add_argument('bamfiles', metavar='bamfile', nargs='+')
    p.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout)
    a = p.parse_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    # Find references
    reference = a.reference
    ref = None

    with a.output as fp:

        bams = (pysam.Samfile(bamfile) for bamfile in a.bamfiles)
        for bam in bams:
            with contextlib.closing(bam):
                if not a.begin or not a.end:
                    positions = [i.pos for i in bam.pileup(reference, max_depth=a.min_reads+1)
                                if i.n >= a.min_reads]
                    if a.begin is None:
                        a.begin = min(positions)
                    if a.end is None:
                        a.end = max(positions) + 1
                    logging.info('Range: %d-%d', a.begin, a.end)
                if not ref:
                    ref = a.reference_fasta.fetch(reference, a.begin, a.end)
                    fp.write('>' + reference + '\n')
                    fp.write(ref + '\n')
                assert(ref)


                for read in bam.fetch(reference, a.begin, a.end):
                    res = []
                    if read.pos > a.begin:
                        res.append('-' * (read.pos - a.begin))
                    ap = read.aligned_pairs
                    # Strip insertions
                    ap = [(q, r) for q, r in ap if r is not None]
                    for q, r in ap:
                        if r < a.begin or r >= a.end:
                            continue
                        if q is None:
                            res.append('-')
                        else:
                            res.append(read.query[q])
                        if r > a.end:
                            break

                    if read.aend < a.end:
                        res.append('-' * (a.end - read.aend))

                    assert(len(''.join(res)) == len(ref))

                    fp.write('>' + read.qname)
                    fp.write('\n' + ''.join(res) + '\n')

if __name__ == '__main__':
    main()
