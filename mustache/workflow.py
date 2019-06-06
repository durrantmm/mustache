from os.path import dirname, join
from snakemake import snakemake


def _workflow(workdir, snakefile, config, cores, memory, unlock):

    snakemake(snakefile, configfile=config, config={'wd': workdir, 'memory': memory}, cores=cores, unlock=unlock)