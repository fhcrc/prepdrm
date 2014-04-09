#!/usr/bin/env python
import argparse
import contextlib
import collections
import csv
import itertools
import logging
import os.path
import operator
import shutil
import subprocess
import sys
import tempfile

from concurrent import futures
import pysam

log = logging.getLogger('prep_drm.contam')

@contextlib.contextmanager
def tempdir(**kwargs):
    td = tempfile.mkdtemp(**kwargs)
    def temppath(*args):
        return os.path.join(td, *args)
    try:
        yield temppath
    finally:
        shutil.rmtree(td)

@contextlib.contextmanager
def blast_database(fasta_file):
    with tempdir(prefix='blast-db') as td:
        cmd = ['makeblastdb',
               '-dbtype', 'nucl',
               '-in', fasta_file,
               '-out', td('blastdb')]
        log.info('Generating BLAST database: %s', ' '.join(cmd))
        subprocess.check_call(cmd)
        yield td('blastdb')

_Sequence = collections.namedtuple('_Sequence',
        ['id', 'seq', 'reference', 'read_group', 'nm', 'alen'])

class Sequence(_Sequence):
    def __init__(self, *args, **kwargs):
        super(Sequence, self).__init__(*args, **kwargs)

    def __len__(self):
        return len(self.seq)

    @property
    def pct_id(self):
        return 1.0 - float(self.nm) / self.alen

def parse_bam(bam_file):
    for read in bam_file:
        rg = read.opt('RG')
        nm = read.opt('NM')
        reference = bam_file.getrname(read.tid)
        yield Sequence(id=read.qname, seq=read.query, reference=reference,
                       read_group=rg, nm=nm, alen=read.alen)

#BlastResult = collections.namedtuple('BlastResult', ['query', 'reference',
    #'pct_id', 'alen', 'nm', 'gapopen_count', 'qstart', 'qend', 'rstart',
    #'rend', 'evalue', 'bitscore'])

BlastResult = collections.namedtuple('BlastResult', ['query', 'reference', 'qlen', 'n_identical', 'qseq', 'sseq'])

def parse_blast6(fp):
    r = csv.reader(fp, delimiter='\t')
    #classes = [((3, 4, 5, 6, 7, 8, 9), int), ((2, 10, 11), float)]
    classes = [((2, 3), int)]
    classes = [(col, cls) for cols, cls in classes for col in cols]

    for row in r:
        for col, cls in classes:
            row[col] = cls(row[col])
        yield BlastResult._make(row)

def blast_reads(sequences, database, max_target_seqs=1):
    count = 0
    with tempfile.NamedTemporaryFile(prefix='blast-query-', suffix='.fasta') as tf:
        for sequence in sequences:
            count += 1
            tf.write('>{0}\n{1}\n'.format(sequence.id, sequence.seq))
        tf.flush()
        log.debug('Wrote %d query sequences to %s', count, tf.name)

        cmd = ['blastn', '-db', database, '-query', tf.name,
               '-outfmt', '6 qseqid sseqid qlen nident qseq sseq',
               '-evalue', '1e-2',
               '-max_target_seqs', str(max_target_seqs)]

        log.debug("BLASTing: %s", ' '.join(cmd))
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

        for record in parse_blast6(p.stdout):
            yield record
        returncode = p.wait()
        if returncode:
            raise subprocess.CalledProcessError(returncode, cmd)

ContaminationRecord = collections.namedtuple('ContaminationRecord',
        ['query', 'read_group', 'expected_reference', 'reference_pct_id', 'top_hit', 'top_hit_pct_id', 'aligned_query', 'aligned_hit'])

def reconcile_reads(blast_result_iter, all_reads):
    seen_reads = set()
    read_dict = {r.id: r for r in all_reads}

    by_query = itertools.groupby(blast_result_iter, operator.itemgetter(0))

    for query_name, results in by_query:
        seen_reads.add(query_name)
        query_info = read_dict[query_name]
        expected_reference = query_info.reference

        top_hit = next(results)
        if (top_hit.reference == expected_reference or
                top_hit.reference.startswith('MB2059_pol') and expected_reference.startswith('MB2059_pol')):
            continue # Unlikely contamination
        else:
            yield ContaminationRecord(query_name,
                    read_group=query_info.read_group,
                    expected_reference=expected_reference,
                    reference_pct_id=query_info.pct_id,
                    top_hit=top_hit.reference,
                    top_hit_pct_id=float(top_hit.n_identical) / top_hit.qlen,
                    aligned_query=top_hit.qseq, aligned_hit=top_hit.sseq)


    # Return anything with *no* hits
    for query_name in set(read_dict) - seen_reads:
        query_info = read_dict[query_name]
        expected_reference = query_info.reference
        yield ContaminationRecord(query_name, read_group=query_info.read_group, expected_reference=expected_reference,
                reference_pct_id=query_info.pct_id, top_hit=None, top_hit_pct_id=None)

def bamfile_task(bam_path, blast_db):
    bam = pysam.Samfile(bam_path)
    with contextlib.closing(bam):
        reads = list(parse_bam(bam))

    log.debug('%s: %d query sequences', bam_path, len(reads))
    blast_result = blast_reads(reads, blast_db)
    contam_result = list(reconcile_reads(blast_result, reads))
    if contam_result:
        log.warn('%s: %d possible contaminants', bam_path, len(contam_result))
    return contam_result

def main():
    p = argparse.ArgumentParser()
    p.add_argument('reference_fasta')
    p.add_argument('bam_paths', nargs='+', metavar='bam_file',
            help="""BAM files to search for contaminants""")
    p.add_argument('-t', '--threads', type=int, default=12)
    p.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout)
    p.add_argument('-v', '--verbose', action='store_const',
            default=logging.INFO, const=logging.DEBUG, dest='level')
    p.add_argument('-q', '--quiet', action='store_const',
            const=logging.WARN, dest='level')
    a = p.parse_args()

    logging.basicConfig(level=a.level,
                        format='[%(levelname)s] %(name)s: %(message)s')

    with blast_database(a.reference_fasta) as db, \
         futures.ThreadPoolExecutor(a.threads) as executor, \
         a.output as ofp:

        w = csv.writer(ofp, lineterminator='\n')
        w.writerow(ContaminationRecord._fields)
        futs = {}
        for bam_path in a.bam_paths:
            f = executor.submit(bamfile_task, bam_path, db)
            futs[f] = bam_path

        while futs:
            done, pending = futures.wait(futs, 1, futures.FIRST_COMPLETED)
            for f in done:
                bam = futs.pop(f)
                log.info('finished %s', bam)
                r = f.result()
                w.writerows(r)

            assert len(futs) == len(pending)
            log.info('%d remaining.', len(pending))

if __name__ == '__main__':
    main()
