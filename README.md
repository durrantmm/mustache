# Installing Mustache

First, install miniconda3. This is an environment management system that should keep everything organized.

Once installed, clone this github directory to some location where it can be stored permanently.

    git clone url
    
Then enter the newly downlaoded mustache directory

    cd mustache
    
And create the new conda environment with

    conda env create -f envs/environment.yml

Now activate the environment with
    
    source activate mustache
    
This is a step that must be repeated whenever using mustache from within this environment.

Now install in an editable form with the command

    pip install --editable .
    
Once complete, you can check to see if mustache installed properly by simply typing

    mustache
    

# Running `mustache align`
Mustache uses BWA to align reads to a reference genome fasta, but it has to properly format the alignment files in order
for them to be compatible with mustache.

We've built in the BWA alignment into the mustache tool, so that you don't have to worry about properly formatting the BAM file.

The format of the argument is

    mustache align (paired|single) FASTQ1 [<FASTQ2>] GENOME OUT_BAM
    
Where the user specifies either paired-end or single-end whole genome sequencing, the 1-2 corresponding FASTQ files,
the reference genome to align to, and the path to the output BAM file.

# Running `mustache find`
The `mustache find` command takes the BAM file produced above and identifies insertion sequences, their genomic loci, 
and the number of reads from the sample supporting the observation.
 
The command is run as

    mustache find (paired|single) [OPTIONS] BAM_FILE OUTPUT_PREFIX

The default options are designed for genomes with 10x coverage, they may need to be adjusted for samples with significantly
higher or lower depth.

The final results will be found in the at `OUTPUT_PREFIX/OUTPUT_PREFIX.stats.tsv`, along with several BAM files with the
relevant aligned reads.