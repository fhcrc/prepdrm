#!/usr/bin/env python
"""
Choose a single reference sequence for each subject.

Sequences are named::
   {SUBJECT_ID}_more...

This script chooses the longest ungapped reference for each subject id.
"""
import argparse
import itertools
import logging
import operator
import subprocess
import sys

from Bio import SeqIO

def consensus(name, sequences):
    sequences = list(sequences)
    if len(sequences) == 1:
        sequence = sequences[0]
        sequence.id = name
        sequence.description = name
        return sequence

    cmd = ['consambig', '-filter', '-name', name]
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE)

    with p.stdin as fp:
        SeqIO.write(sequences, fp, 'fasta')
    with p.stdout as fp:
        try:
            return SeqIO.read(fp, 'fasta')
        finally:
            return_code = p.wait()
            if return_code:
                raise subprocess.CalledProcessError(return_code, cmd)

def ungapped_len(sequence):
    return sum(1 for c in str(sequence.seq) if c != '-')

def choose_longest(name, sequences):
    sequences = list(sequences)
    #if name == '5111318':
        #logging.info('%s: %s %s', name, [i.id for i in sequences],
                     #[ungapped_len(i) for i in sequences])
    sequence = max(sequences, key=ungapped_len)
    old_id = sequence.id
    sequence.id = name
    sequence.description = old_id
    return sequence

def split_id(sequence):
    return sequence.id.split('_', 1)[0]

def main():
    p = argparse.ArgumentParser(description=__doc__)

    p.add_argument('-i', '--input', type=argparse.FileType('r'),
            default=sys.stdin, help="""Input FASTA file [default: stdin]""")
    p.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help="""Input FASTA file [default: stdout]""")

    a = p.parse_args()
    logging.basicConfig(level=logging.INFO)

    with a.input as fp:
        reads = sorted(SeqIO.parse(fp, 'fasta'), key=operator.attrgetter('id'))

    grouped = itertools.groupby(reads, split_id)
    #consensus_sequences = (consensus(g, seqs) for g, seqs in grouped)
    longest = (choose_longest(g, seqs) for g, seqs in grouped)
    SeqIO.write(longest, a.output, 'fasta')


if __name__ == '__main__':
    main()
