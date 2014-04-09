from setuptools import setup, find_packages
from Cython.Build import cythonize

ext_modules = cythonize('prep_drm/*.pyx')

setup(
    name='prep_drm',
    author='Connor McCoy',
    author_email='cmccoy@fhcrc.org',
    packages=find_packages(),
    #package_data={'prep_drm': ['imgt/data/*.fasta']},
    entry_points={
        'console_scripts': [
            'bam_consensus.py = prep_drm.scripts.bam_consensus:main',
            'bam_to_msa.py = prep_drm.scripts.bam_to_msa:main',
            'classify_mutations.py = prep_drm.scripts.classify_mutations:main',
            'errors_in_region.py = prep_drm.scripts.errors_in_region:main',
            'hypermutation.py = prep_drm.scripts.hypermutation:main',
            'k65r_poisson_binomial.py = prep_drm.scripts.k65r_poisson_binomial:main',
            'make_sample_json.py = prep_drm.scripts.make_sample_json:main',
            'recode_homopolymer_regions.py = prep_drm.scripts.recode_homopolymer_regions:main',
            'stats_by_readgroup.py = prep_drm.scripts.stats_by_readgroup:main',
            'translate_gff3.py = prep_drm.scripts.translate_gff3:main',
    ]},
    test_suite='prep_drm.test.suite',
    ext_modules=ext_modules,
    install_requires=[])
