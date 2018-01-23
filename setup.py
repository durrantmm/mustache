from setuptools import setup, find_packages

setup(
    name="mustache",
    version='0.1',
    description='Tool to identify insertion sequences from whole-genome sequencing data.',
    url='https://github.com/durrantmm/mustache',
    author="Matt Durrant",
    author_email="mdurrant@stanford.edu",
    license="MIT",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'click',
        'pysam',
        'snakemake',
        'biopython',
        'scipy',
        'pandas',
        'statsmodels',
        'numpy'
    ],
    zip_safe=False,
    entry_points = {
        'console_scripts': [
            'mustache = mustache.main:cli'
            # more script entry points ...
        ],
}
)
