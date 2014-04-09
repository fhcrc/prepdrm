#!/usr/bin/env python
import argparse
import collections
import csv
import json
import os.path
import sys

def main():
    p = argparse.ArgumentParser()
    p.add_argument('plate_json', type=argparse.FileType('r'))
    p.add_argument('sample_meta', type=argparse.FileType('r'))
    p.add_argument('-o', '--outfile', type=argparse.FileType('w'),
            default=sys.stdout)
    a = p.parse_args()

    with a.plate_json as fp:
        plate_info = json.load(fp, object_pairs_hook=collections.OrderedDict)

    with a.sample_meta as fp:
        # Drop comments
        lines = (i for i in fp if not i.startswith('#'))
        sample_meta = list(csv.DictReader(lines))

    result = plate_info.copy()
    result['samples'] = []

    for row in sample_meta:
        sample_sam_tags = result['sam_tags'].copy()
        sample_sam_tags['ID'] = row['sample_id']
        sample_sam_tags['PU'] = row['barcode']
        sample_sam_tags['SM'] = row['subject']

        sample = collections.OrderedDict()
        sample['fastq_path'] = os.path.join(result['split_fastq_directory'],
                                            row['sample_id'] + '.fastq.bz2')
        if not os.path.exists(sample['fastq_path']):
            print >> sys.stderr, "DNE:", sample['fastq_path']
        sample['plate'] = int(row['plate'])
        sample['sample_id'] = row['sample_id']
        sample['is_control'] = row['is_control'].upper() == 'TRUE'
        sample['direction'] = row['direction']
        sample['barcode_id'] = row['barcode_id']
        sample['reference_id'] = row['reference_id']
        sample['sam_tags'] = sample_sam_tags

        result['samples'].append(sample)


    with a.outfile as ofp:
        json.dump(result, ofp, indent=2)

if __name__ == '__main__':
    main()
