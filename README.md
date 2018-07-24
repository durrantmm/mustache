# Installing Mustache

First, install miniconda3. This is an environment management system that should keep everything organized.

Once installed, clone this github directory to some location where it can be stored permanently.

    git clone https://github.com/durrantmm/mustache.git
    
Then enter the newly downlaoded mustache directory

    cd mustache
    
Install the mustache conda environment and other dependencies with

    bash install.sh

Now activate the environment with
    
    source activate mustache
    
This is a step that must be repeated whenever using mustache from within this environment.

Now install mustache from with the command

    pip install --editable .
    
Once complete, you can check to see if mustache installed properly by simply typing

    mustache
   
This can then be called from anywhere on the file system while in the `mustache` conda environment.

# Running `mustache alignbwa`
Mustache uses BWA to align reads to a reference genome fasta, but it has to properly format the alignment files in order
for them to be compatible with certain mustache tools, such as alignment_anchord assembly.

We've built in the BWA alignment into the mustache tool, so that you don't have to worry about properly formatting the BAM file.

The format of the argument is

    python alignbwa.py GENOME OUT_BAM FASTQ1 [<FASTQ2>] 
or
  
    mustache align GENOME OUT_BAM FASTQ1 [<FASTQ2>] 
    
Whether or not the alignment is paired-end is determined by whether one or two fastq files are provided.

# Running `mustache find`
The `mustache find` command takes the BAM file produced above and identifies insertion sequences, their genomic loci, 
and the number of reads from the sample supporting the observation.
 
The command is run as

    mustache find (paired|single) [OPTIONS] BAM_FILE GENOME OUTPUT_PREFIX

The default options are designed for genomes with 10x coverage, they may need to be adjusted for samples with significantly
higher or lower depth. 

## `mustache find` output
The final results will be found in the specified output directory.

    *.insertions_seqs.tsv - Detailed table of insertion sequences identified.
    *.insertions_seqs.fasta - A fasta file of all the identified insertion sequences.
