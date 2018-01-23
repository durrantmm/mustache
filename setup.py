from setuptools import setup

setup(
    name="mustache",
    version='0.1',
    description='Tool to identify insertion sequences from whole-genome sequencing data.',
    url='https://github.com/durrantmm/mustache',
    py_modules=['hello'],
    install_requires=[
        'click',
        'pysam',
        'snakemake',
        'biopython',
        'scipy',
        'pandas',
        'statsmodels',
        'numpy',
    ],
    entry_points='''
        [console_scripts]
        mustache=mustache.main:cli
    ''',
)
