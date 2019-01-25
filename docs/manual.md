[Back to main page](../README.md)  

# User Manual 
This user manual is designed as a reference for users who have questions about *mustache*. It includes details about
how to execute important commands, the parameters available, brief descriptions of how the command works, and 
descriptions of the output files.

## Primary *mustache* commands
We now describe each of the primary *mustache* commands. These commands are at the center of the *mustache* workflow,
and should not be ignored when doing a complete analysis.

### `findflanks`
![alt text](img/findflanks.png)

This command finds insertion sites and reconstructs the flanks of inserted sequence.

#### `findflanks`: Input and parameters
The `findflanks` command takes a BAM file as input. A local alignment algorithm, such as 
[BWA MEM](http://bio-bwa.sourceforge.net/) should be used to generate this BAM file.

You can run the command as

    mustache findflanks BAMFILE
    

Additional parameters include:

    --min_softclip_length, -minlen
    --min_softclip_count, -mincount
    --min_count_consensus, -mcc
    --min_alignment_quality, -minq
    --min_softclip_ratio, -minratio
    --min_alignment_inner_length, -minial
    --min_distance_to_mate, -mindist
    --max_indel_ratio, -maxir

The parameter `min_softclip_length` takes an integer. 
For a softclipped site to be considered, there must be at least one softclipped read of this length. 
The default used is 8.
 
The parameter `min_softclip_count` takes an integer.
It requires that for a clipped site to be considered it must be supported by at least this number of reads.
The default used is 4.

The parameter `min_count_consensus` takes an integer.
When generating a consensus sequence for the insertion flank, it requires that each site must have at least this many 
reads supporting it.
The default used is 2.

The parameter `min_alignment_quality` takes an integer.
This determines the minimum mapping quality of each alignment in order for it to be considered.
The default used is 20.

The parameter `min_softclip_ratio` takes a float between 0 and 1.0.
It sets a minimum for the ratio of clipped reads to all reads at a given site.
This filters out very small indels that could disrupt the flank-pairing step.
The default used is 0.15.

The parameter `min_alignment_inner_length` takes an integer.
It filters out reads that are softclipped on both ends, and where the inner aligned portion is shorter than this length.
Filtering out such reads removes a potential false positives.
The default is 21, this parameter should not be changed.

The parameter `min_distance_to_mate` takes an integer.
It filters out clipped-sites that have no other candidate flank pair nearby. 
In the `pairflanks` step, we only pair flanks with each other if they are within 20 bp of 
each other. Filtering out sites without a candidate pair can dramatically speed up the 
analysis. However, if you want to keep all candidate insertion sites, you can increase this
number to some very high integer. The default is 22.

#### `findflanks`: Description of implementation
The `findflanks` algorithm works by identifying candidate insertion sites by searching for clipped-end sites in locally 
aligned reads. To generate a consensus sequence of the candidate flank, we use a trie-based approach intended to filter 
out spurious reads when building a consensus sequence.  This also allows our algorithm to distinguish between multiple 
insertions at a single site, which may be observed in a metagenomic sample, for example (analysis of metagenomic samples
is in development).

Briefly, we build a sequence Trie of all of the reads found clipped at a given site. We then generate all of the 
unique, complete sequences by traversing the trie. We then cluster these sequences together by sequence similarity, 
and process the different clusters separately.

We then generate a consensus sequence for a given insertion flank by traversing down the trie, choosing the base with 
the most support at each step. Evidence for a given base is calculated as the sum of all base-pair quality scores 
originally reported in the FASTQ file.

The consensus sequence ends when the number of reads supporting a position drops below `--min_count_consensus`.

#### `findflanks`: Output file format
By default, `findflanks` will write to a file named `mustache.findflanks.tsv` (The name can be changed with the 
`--output_file` parameter). The columns of this file are described as follows:

1. `flank_id` - An arbitrary identifier used internally by *mustache*.
2. `contig` - The name of the contig where the flank was identified.
3. `pos` - The exact base pair position of where the clipped-ends begin, the presumed start of the insertion. 
This location is 0-based with respect to the reference genome used.
4. `orient` - The orientation of the flank. `5p` refers to the 5' end of an insertion, which suggests this flank runs 
from the 5' to the 3' direction with respect to the reference genome. `3p` refers to the 3' end of an insertion, which 
suggests the flank runs from the 3' to the 5' direction with respect to the reference genome.
5. `softclip_count_5p` - The total number of reads clipped from the 5' to 3' direction at this site.
6. `softclip_count_3p` - The total number of reads clipped from the 3' to 5' direction at this site.
7. `runthrough_count` - The total number of reads that run through the insertion site without being clipped.
8. `small_insertion_count_5p` - The total number of reads that include a small insertion from the 5' to 3' direction at
this site. If this number is high, it indicates that the detected insertion is most likely a small indel. `findflanks` 
uses this to filter out false positives.
9. `small_insertion_count_3p` - Same as above, but for small insertions running from the 3' to the 5' direction.
10. `deletion_count` - The number of reads with a deletion at the insertion site, also used to filter out false 
positives.
11. `upstream_deletion_count` - The number of reads that contain a deletion one base pair upstream of this insertion 
site. Also
useful for removing false positives
12. `downstream_deletion_count` - Same as above, but for reads with deletions one base pair downstream of this insertion
site.
13. `total_count` - The sum of the `softclip_count_5p`, `softclip_count_3p`, `runthrough_count`, 
`small_insertion_count_5p`, `small_insertion_count_3p`, and `deletion_count` columns.
14. `consensus_softclip_count` - The number of clipped reads at this site that specifically support the 
consensus flank sequence on this line. This will be used as the `softclip_count_5p` and `softclip_count_3p` columns
in the file generated by running the `pairflanks` command.
15. `consensus_seq` - The consensus sequence for the insertion flank identified at this site.
    
This intermediate file is not particularly valuable to users on its own, and it is primarily intended to be used
internally by *mustache* in the next step.

### `pairflanks`
This command pairs insertion flanks with each other to represent the 5' and 3' flanks of a candidate insertion.

#### `pairflanks`: Input and parameters
You can run the `pairflanks` command as

    mustache pairflanks FLANKSFILE BAMFILE REFERENCE_GENOME

Where `FLANKSFILE` is the output of the `findflanks` command, `BAMFILE` is the binary sequence alignment file used as 
input for `findflanks`, and `REFERENCE_GENOME` is the reference genome that `BAMFILE` was aligned to.

Additional parameters include:

    --max_direct_repeat_length, -maxdr
    --min_alignment_quality, -minq
    --min_alignment_inner_length, -minial
    --max_junction_spanning_prop, -maxjsp

The parameter `max_direct_repeat_length` takes an integer. This species the maximum distance that oppositely-oriented
insertion flanks can be from each other in order to consider pairing them together. Since insertions often cause direct
repeats at the insertion site, the location of the insertion flanks are often separated by several base pairs. 
The default for this parameter is 20. This means that our algorithm by default will NOT identify true insertions that 
create direct repeats that exceed 20 base pairs.

The parameter `min_alignment_quality` takes an integer. This determines the minimum mapping quality of each alignment
in order for it to be considered. The default used is 20. This should be the same as the `--min_alignment_quality`
parameter in `findflanks`.

The parameter `min_alignment_inner_length` takes an integer.
It filters out reads that are softclipped on both ends, and where the inner aligned portion is shorter than this length.
Filtering out such reads removes a potential false positives.
The default is 21, this parameter should not be changed, and it should be the same as the value of 
`--min_alignment_inner_length` used in `findflanks.`

The parameter `max_junction_spanning_prop` takes an integer. This is used to filter out low-confidence insertion pairs,
or insertions that are occurring in duplicated regions. We expect that very few reads will fully span the insertion
junction, as that would indicate that the insertion site does not contain the insertion somewhere in the isolate. 
If the number of reads spanning the insertion junction without being clipped exceeds this proportion of the total reads
at the site, then it will be ignored. By default, this parameter is 0.15.

#### `pairflanks`: Description of implementation
The `pairflanks` command uses a variety of techniques to pair flanks with each other. It first filters insertion flanks 
by the `--max_direct_repeat_length` parameter. It then does pairwise comparisons between all nearby, oppositely-oriented 
insertion flanks to determine if they share inverted repeats at their termini, a common feature of many prokaryotic 
insertion sequences. It then prioritizes all pairs by the following parameters, in order: 
1. The length of the inverted repeat identified, if any. Pairs with longer shared inverted repeats receive higher 
priority.
2. The difference in the number of softclipped reads that support the insertion flanks. 
Pairs that have a similar number of reads supporting each flank are given higher priority.
3. The difference in the length of the recovered insertion flanks. 
If the two insertion flanks have very similar lengths, they are given higher priority than pairs with disparate flank 
lengths.

It then filters sites with too many junction-spanning reads (see `max_junction_spanning_prop`), and infers the identity
of the direct repeat created by the insertion, if any. We realize that these are somewhat arbitrary filters, and we 
welcome suggestions from the users on how this might be improved.

#### `pairflanks`: Output file format
By default, `pairflanks` will write to a file named `mustache.pairflanks.tsv` (The name can be changed with the 
`--output_file` parameter). The columns of this file are described as follows:

1. `pair_id` - An arbitrary identifier used to identify the candidate insertion. This is an important identifier that
is used to relate the output of the `pairflanks` command to the output of the various `inferseq` commands.
2. `contig` - The name of the contig where the pair was identified.
3. `pos_5p` - The exact base pair position of where the clipped-ends of the 5' insertion flank begin, the presumed 
start of the insertion. This location is 0-based with respect to the reference genome used.
3. `pos_3p` - The exact base pair position of where the clipped-ends of the 3' insertion flank begin, the presumed 
end of the insertion. This location is 0-based with respect to the reference genome used.
4. `softclip_count_5p` - The number of clipped reads at this site that support the 5' insertion flank. This is taken 
directly from the `consensus_softclip_count` of the `findflanks` output file.
5. `softclip_count_4p` - The number of clipped reads at this site that support the 3' insertion flank. This is also 
taken directly from the `consensus_softclip_count` of the `findflanks` output file.
6.  `total_count_5p` - The total number of reads found at the insertion site for the 5' flank.
7. `total_count_3p` - The total number of reads found at the inseriton site for the 3' flank.
8. `spanning_count` - The total reads that are found to span the 5' and 3' insertion sites by 10 bp in both directions.
9. `has_IR` - A boolean True/False indicating whether or not the flank pair contains inverted repeats at their end, a 
common feature of many inserted elements.
10. `IR_length` - The length of the detected inverted repeat.
11. `IR_5p` - The inverted repeat in the 5' flank.
12. `IR_3p` - The inverted repeat in the 3' flank.
13. `seq_5p` - The consensus sequence of the 5' flank.
14. `seq_3p` - The consensus sequence of the 3' flank.
15. `direct_repeat_reference` - The sequence of the direct repeat created, using the reference genome to infer this
direct repeat. 
16. `direct_repeat_reads_consensus` - The sequence of the direct repeat created, but using information from the reads 
themselves to determine the a more accurate direct repeat. This can be useful if the direct repeat of the isolate
actually contains some base pair substitutions that do not exist in the reference genome.

This file will likely have intrinsic value to the user, as they indicate candidate insertions. However, these candidates
must be further investigated before any final conclusions about their identity can be made. In some cases, for example,
these may actually represent genomic inversions, not insertions. To infer the identity of these insertions, several
inference approaches are implemented and described in the following sections.

### `inferseq-reference`
![alt text](img/inferseqreference.png)

This command infers the identity of insertions by aligning the flanks of candidate pairs to a reference genome.

#### `inferseq-reference`: Input and parameters
The `inferseq-reference` command takes a the output of the `pairsfile` command and a reference genome as input.

While *mustache* will automatically index the reference genome if it is not already indexed, we recommend that you index 
the reference genome beforehand with the command:

    bowtie2-build -o 0 -q INFERSEQ_REFERENCE INFERSEQ_REFERENCE

Where INFERSEQ_REFERENCE is the reference genome of interest.

[BWA MEM](http://bio-bwa.sourceforge.net/) should be used to generate this BAM file.

You can run then run the command as

    mustache inferseq-reference PAIRFLANKS INFERSEQ_REFERENCE
    
Additional parameters include:

    --min_perc_identity, -minident
    --max_internal_softclip_prop, -maxclip
    --max_inferseq_size, -maxsize
    --min_inferseq_size, -minsize
    --keep-intermediate/--no-keep-intermediate


The parameter `min_perc_identity` takes an float between 0 and 1. When aligning candidate insertion flanks to a the
reference genome, it will only consider alignments that exceed this percentage identity with the reference. The default
for this parameter is 0.95. 

The parameter `max_internal_softclip_prop` takes an float between 0 and 1. This is an additional candidate insertion
flank alignment filter. If the aligned flanks are internally clipped by a proportion of their total length that exceeds
this number, then the alignment is excluded. For example, imagine the alignment:

    -------->          <----------
    flank1                  flank2

Imagine the vertical bar `|` indicating that one of the alignments is clipped at a specific site, such as:
    
    -------->         <---|-------
    flank1                  flank2

If the proportion of the clipped end of the flank2 alignment exceeds 0.05 (by default), then this alignment will be
excluded.

The parameter `max_inferseq_size` excludes inferred sequences that are larger than `max_inferseq_size`. The default for
this parameter is 500 kilobase pairs, but we urge users to be cautious when working with very large insertions, as these
may be false positives.

The parameter `min_inferseq_size` excludes inferred sequences that are smaller than `min_inferseq_size`. The default for
this parameter is 1 base pair, but in reality *mustache* is not well suited to identify insertions that are smaller than
the read length of the sequencing library.

The `--keep-intermediate/--no-keep-intermediate` will keep determine whether or not the intermediate alignment file will
be deleted or kept, which can be useful for debugging purposes.

#### `inferseq-reference`: Description of implementation
#### `inferseq-reference`: Output file format

### `inferseq-assembly`
![alt text](img/inferseqassembly.png)

This command finds insertion sites and reconstructs the flanks of inserted sequence.

#### `inferseq-assembly`: Input and parameters
#### `inferseq-assembly`: Description of implementation
#### `inferseq-assembly`: Output file format


### `inferseq-overlap`

#### `inferseq-overlap`: Input and parameters
#### `inferseq-overlap`: Description of implementation
#### `inferseq-overlap`: Output file format

### `inferseq-database`
![alt text](img/inferseqdatabase.png)

#### `inferseq-database`: Input and parameters
#### `inferseq-database`: Description of implementation
#### `inferseq-database`: Output file format



## Additional/ Experimental *mustache* Commands
### `formatbam`

#### `formatbam`: Input and parameters
#### `formatbam`: Description of implementation
#### `formatbam`: Output file format

### `recall`

#### `recall`: Input and parameters
#### `recall`: Description of implementation
#### `recall`: Output file format

### `extendpairs`

#### `extendpairs`: Input and parameters
#### `extendpairs`: Description of implementation
#### `extendpairs`: Output file format

[Back to main page](../README.md)