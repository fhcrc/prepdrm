from datetime import datetime
import glob
import itertools
import json
import os
import os.path

#from bioscons.slurm import SlurmEnvironment
from bioscons.fast import fast

from SCons.Script import Builder, Environment


env = Environment()
env.AppendENVPath('PATH', '/app/bin')
env.AppendENVPath('PATH', '/home/matsengrp/local/encap/pandoc-1.9.4.2/bin')
env.PrependENVPath('PATH', '/home/matsengrp/local/encap/R-3.0.2/bin')
env.PrependENVPath('PATH', './bin')
env.PrependENVPath('CLASSPATH', 'lib/picard-1.84.jar:lib/sam-1.84.jar')
env['JVM_OPTS'] = '-XX:MaxPermSize=256m -Xmx512m'
env['outdir'] = 'output'
#env['ENV']['SLURM_ACCOUNT'] = 'overbaugh_j'
env.VariantDir('output', src_dir='.')
fast(env)

def stripext(s, basename=False):
    r = os.path.splitext(s)[0]
    if basename:
        return os.path.basename(r)
    return r


# Builders
############
sam_to_bam = Builder(action='samtools view -S -b -o $TARGET $SOURCE',
                     suffix='.unsorted.bam', src_suffix='.sam')
env['BUILDERS']['SamToBam'] = sam_to_bam

def picard_merge(source, target, env, for_signature):
    sources = map(str, source)
    source_str = ' '.join('I=' + i for i in sources)

    cmd = ('java $JVM_OPTS net.sf.picard.sam.MergeSamFiles '
           'USE_THREADING=true MERGE_SEQUENCE_DICTIONARIES=true '
           'VALIDATION_STRINGENCY=LENIENT ' +
           source_str + ' O=' + str(target[0]))
    return cmd
env['BUILDERS']['MergeBam'] = Builder(generator=picard_merge)
index_bam = Builder(action='samtools index $SOURCE',
        suffix='.bam.bai', src_suffix='.bam')
env['BUILDERS']['IndexBam'] = index_bam
index_fasta = Builder(action='samtools faidx $SOURCE',
        suffix='.fasta.fai', src_suffix='.fasta')
env['BUILDERS']['FaIndex'] = index_fasta
bam_to_msa = Builder(action='bam_to_msa.py -r $REFERENCE $REF_FASTA $SOURCE -o $TARGET',
        src_suffix='.bam', suffix='.msa.fasta')
env['BUILDERS']['BamToMsa'] = bam_to_msa

# Archiving
def archive_analysis(env, source):
    archive_dir = 'prep_drm_analysis_{0}'.format(datetime.now().strftime('%Y%m%d'))
    archive_dir = env.get('archive_dir', archive_dir)

    sources = sorted([str(i) for i in source])

    to_zip = []
    for orig_dir, files in itertools.groupby(sources, os.path.dirname):
        assert orig_dir.startswith('output'), orig_dir
        dest_dir = os.path.join(archive_dir, (orig_dir + '/')[7:])
        to_zip.extend(env.Install(dest_dir, list(files)))

    z = env.Zip(archive_dir + '.zip', to_zip)
    # Delete intermediates
    deletion = env.Command('.nonexistant-file', to_zip, 'rm -rf {0}'.format(archive_dir))
    env.Depends(deletion, z)
    return z, deletion
env.AddMethod(archive_analysis, 'ArchiveAnalysis')
############

# Start Building
################



# Build references
env['REF_FASTA_ALN'] = 'references/prep_refs_MB2059.fasta'
env['MUTATION_DIRECTIONS'] = 'metadata/location_directions.csv'
ref_fasta, = env.Command('$outdir/prep_references.fasta',
        '$REF_FASTA_ALN',
        'seqmagick convert --ungap $SOURCE $TARGET')
env['REF_BOUNDS'] = 'references/prep_refs_MB2059.bounds.json'
ref_fasta_index = env.FaIndex(ref_fasta)
env['REF_FASTA'] = ref_fasta
# Contamination set
contam_screen = env.Command('$outdir/contam_screen_references.fasta',
        [ref_fasta, 'references/HXB2_pol_amplicon.fasta', 'references/overbaugh_refs.fasta'],
        'cat $SOURCES | seqmagick convert --ungap - $TARGET')

merged_metadata, = env.Command('$outdir/prep_merged_metadata.csv',
        env.Glob('metadata/plate[0-9].csv'),
        'stack_csvs.r $SOURCES > $TARGET')

################
# Data Loading #
################
def load_samples():
    plates = []
    for f in glob.iglob('metadata/*_metadata.json'):
        with open(f) as fp:
            plates.append(json.load(fp))
    return plates


def generate_references(plates):
    with open(env['REF_BOUNDS']) as fp:
        bounds = json.load(fp)
    reference_ids = set(sample['reference_id'] for plate in plates
                        for sample in plate['samples'])

    references = {}

    for reference_id in reference_ids:
        gff3 = env.Command('$outdir/references/{0}.gff3'.format(reference_id),
                ['$REF_FASTA_ALN', 'references/MB2059_pol.gff3'],
                'translate_gff3.py -q {0} $SOURCES -o $TARGET'.format(reference_id))[0]
        fasta = env.Command('$outdir/references/{0}.fasta'.format(reference_id),
                ['$REF_FASTA_ALN'],
                'seqmagick convert --pattern-include "^{0}$" --ungap $SOURCE $TARGET'.format(reference_id))[0]
        references[reference_id] = {'fasta': fasta, 'gff3': gff3, 'faidx': None,
                                    'bounds': tuple(bounds[reference_id])}

    return references


def analyze_sample(env, sample, references):
    targets = ['$outdir/indiv/{bn}{ext}'.format(bn=sample['sample_id'], ext=ext)
               for ext in ('_nt.fasta', '.unsorted.bam')]
    fastq = sample['fastq_path']
    ref = references[sample['reference_id']]
    sources = [ref['fasta'], fastq]
    start, end = ref['bounds']
    bounds_args = '-b {0} -e {1}'.format(start, end)

    e = env.Clone()
    pairwise_aln, bam = e.Command(targets, sources,
            'codon-sw -hpfs -7 -q $SOURCES --fasta-pairs $TARGETS ' +
            bounds_args)
    env.Precious([pairwise_aln, bam])
    #pairwise_aln, bam = e.Command(targets, sources,
            #'codon-sw -hpfs -7 -q  $SOURCES --fasta-pairs $TARGETS')
    e.Depends(bam, references[sample['reference_id']]['faidx'])

    tags = ' '.join('{0}="{1}"'.format(k, v)
                    for k, v in sample['sam_tags'].iteritems()
                    if k in ('ID', 'DT', 'SM', 'LB', 'PL', 'PU'))
    sort_command = ('java $JVM_OPTS '
         'net.sf.picard.sam.AddOrReplaceReadGroups ' +
         tags +
         ' I=$SOURCE O=$TARGET SORT_ORDER=coordinate '
         ' VALIDATION_STRINGENCY=SILENT QUIET=true CREATE_INDEX=false')

    sorted_bam = e.Command(str(bam).replace('.unsorted', ''),
            bam, sort_command)

    ## Index
    bam_index_file = e.IndexBam(sorted_bam)
    # MSA
    e['REFERENCE'] = sample['reference_id']
    bam_msa, = e.BamToMsa(sorted_bam)
    env.Depends(bam_msa, ref_fasta)
    env.Depends(bam_msa, ref_fasta_index)
    env.Depends(bam_msa, bam_index_file)
    return {'pairwise': pairwise_aln,
            'bam': sorted_bam,
            'msa': bam_msa,
            'index': bam_index_file}


############
# Analysis #
############
plates = load_samples()
references = generate_references(plates)

merged_gff3, = env.Command('$outdir/prep_references.gff3',
        [i['gff3'] for i in references.itervalues()],
        'cat $SOURCES > $TARGET')

pairwise_alignments = []
multiple_alignments = []
sorted_bams = []
to_archive = []
classified = []
classified_largecontext = []

for plate in plates:
    label = plate['label']
    samples = plate['samples']
    e = env.Clone()
    e['outdir'] = os.path.join(env['outdir'], label)
    e['label'] = label
    plate_bams = []
    plate_indexes = []
    for sample in samples:
        r = analyze_sample(e, sample, references)
        plate_bams.append(r['bam'])
        sorted_bams.append(r['bam'])
        pairwise_alignments.append(r['pairwise'])
        multiple_alignments.append(r['msa'])

    pmerged, = e.MergeBam('$outdir/${label}.bam', plate_bams)

    e_blast = e.Clone()
    cpus = 1
    e_blast['BLAST_THREADS'] = cpus
    blast_res, = e_blast.Command('$outdir/${label}_contam_hits.csv',
                            e.Flatten([contam_screen, plate_bams]),
                            'blast_against_refs.py $SOURCES -o $TARGET '
                            '--threads $BLAST_THREADS -q')

    pmerged, = e.Command('$outdir/${label}_prunelowid.bam',
            [pmerged, blast_res],
            'drop_low_identity.py $SOURCE $TARGET --blast-contaminants=${SOURCES[1]}')
    pmerged, = e.Command('$outdir/${label}_prunelowid_recodehomopol.bam',
                        [ref_fasta, merged_gff3, pmerged],
                        'recode_homopolymer_regions.py $SOURCES $TARGET && samtools index $TARGET')

    csv_stats, = e.Command('$outdir/${label}_classified.csv',
            [pmerged, merged_gff3],
            'classify_mutations.py $SOURCES -o $TARGET')
    classified.append(csv_stats)
    csv_stats_largecontext = e.Command('$outdir/${label}_classified_largecontext.csv',
            [pmerged, merged_gff3],
            'classify_mutations.py $SOURCES -o $TARGET --context=9')
    classified_largecontext.append(csv_stats_largecontext)


sorted_bams.sort(key=str)
env.Alias('bams', sorted_bams)

merged, = env.MergeBam('$outdir/prep_merged.bam', sorted_bams)
merged_index, = env.IndexBam(merged)

## Run stats

csv_stats, = env.Command('$outdir/prep_merged_classified.csv',
        classified,
        'csvstack $SOURCES > $TARGET')
csv_stats_largecontext, = env.Command('$outdir/prep_merged_classified_largecontext.csv',
        classified_largecontext,
        'csvstack $SOURCES > $TARGET')
by_sample_largecontext = env.Command(
    '$outdir/prep_merged_classified_largecontext_bysample.csv',
    [csv_stats_largecontext, merged_metadata],
    'sequence_classification_by_sample.r $SOURCES $TARGET')
plots, counts = env.Command(['$outdir/prep_merged_classified.pdf',
                             '$outdir/prep_merged_classified_counts.csv'],
                            csv_stats,
                            'analyze_results.r $SOURCES $TARGETS')

agg_aa_counts, agg_aa_wide = env.Command(['$outdir/prep_merged_aa_counts.csv',
                                          '$outdir/prep_merged_aa_counts_wide.csv'],
                                         [merged_metadata, csv_stats],
                                         'aggregate_aa_per_sample.r $SOURCES $TARGETS')

alen_stats, = env.Command('$outdir/prep_alignment_stats.csv',
        merged,
        'stats_by_readgroup.py $SOURCE -o $TARGET')

alen_stats_plot, = env.Command('$outdir/prep_alignment_stats.pdf',
        [alen_stats, merged_metadata],
        'plot_alignment_stats.r $SOURCES $TARGET')

p_values = env.Command(['$outdir/prep_merged_fisher.csv',
                        '$outdir/prep_merged_fisher_wide.csv',
                        '$outdir/prep_merged_fisher.pdf'],
        [merged_metadata, counts],
        'enrichment_fisher_test.r $SOURCES $TARGETS')

for subset in ('all', 'primary'):
    to_archive.extend(env.Command('$outdir/prep_subject_results_wide_{0}.csv'.format(subset),
        [merged_metadata, p_values[0], agg_aa_wide],
        'prep_subject_results_wide.r {0} $SOURCES $TARGET'.format(subset)))

to_archive.extend([csv_stats, csv_stats_largecontext, plots, counts, p_values,
                   alen_stats_plot, by_sample_largecontext])

env.Default('$outdir')


## Archive and zip
archive_dir = 'prep_drm_analysis_{0}'.format(datetime.now().strftime('%Y%m%d'))
z = env.ArchiveAnalysis(env.Flatten([to_archive]))
env.Alias('zip', z)
