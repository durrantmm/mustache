from setuptools import setup

setup(
    name="mustache",
    version='0.1',
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
        mustache=main:cli
    ''',
)
