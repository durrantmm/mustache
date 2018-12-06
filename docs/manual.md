[Back to main page](../README.md)  

We now describe each of the *mustache* commands.
## `findflanks`
This command finds insertion sites and reconstruct flanks of inserted sequence

The findflanks command takes a BAM file as input. A local alignment algorithm, such as [BWA MEM](http://bio-bwa.sourceforge.net/)
should be used to generate this BAM file.

You can run the command as

    mustache findflanks BAMFILE
    
The findflanks algorithm works by identifying candidate insertion sites by searching for clipped-end sites in locally aligned reads.
It then filters sites according to a variety of user-specified parameters, such as

    -minlen, --min_softclip_length
    -mincount, --min_softclip_count
    -mcc, --min_count_consensus
    -minq, --min_alignment_quality
    -minratio, --min_softclip_ratio
    -minial, --min_alignment_inner_length

The parameter `min_softclip_length` takes an integer. 
For a softclipped site to be considered, there must be at least one softclipped read of this length. 
The default used is 4.
 
The parameter `min_min_softclip_count` takes an integer.
It requires that for a clipped site to be considered it must be supported by at least this number of reads.
The default used is 4.

The parameter `min_count_consensus` takes an integer.
When generating a consensus sequence for the insertion flank, it requires that each site must have at least this many reads supporting it.
The default used is 2.

The parameter `min_alignment_quality` takes an integer.
This determines the minimum mapping quality of each alignment in order for it to be considered.
The default used is 20.

The parameter `min_softclip_ratio` takes a float between 0 and 1.0.
It sets a minimum for the ratio of clipped reads to non-clipped reads at a given site.
This filters out very small indels that could disrupt the flank-pairing step.

The parameter `min_alignment_inner_length` takes an integer.
It filters out reads that are softclipped on both end, and where the inner aligned portion is shorter than this length.
Filtering out such reads removes a lot of false positives.
The default is 21, this parameter should not be changed.

Additionally, we filter insertion flanks to only include sites with nearby, oppositely-oriented insertion flanks.
This reduces computation time considerably.


To generate a consensus sequence of the candidate flank, we use a trie-based approach intended to filter out spurious reads when
building a consensus sequence.  This also allows our algorithm to distinguish between multiple insertions at a single site,
which may be observed in a metagenomic sample, for example.

Briefly, we build a sequence Trie of all of the reads found clipped at a given site. We then generate all of the unique, complete
sequences by traversing the trie. We then cluster these sequences together by sequence similarity, and process the different clusters
separately.

We then generate a consensus sequence for a given insertion flank by traversing down the trie, choosing the base with the most support at each step.
Evidence for a given base is calculated as the sum of all base-pair quality scores originally reported in the FASTQ file.
The consensus sequence ends when the number of reads supporting a position drops below `min_count_consensus`.

Here is an example output file:

    flank_id        contig  pos     orient  softclip_count  runthrough_count        consensus_seq
    1      NC_007795.1     95766   L       23      28      GATAAGTAGAAATGGTAAAAACATTGTATAGCATTTTACACAGGAGTCTGGACTTGACT
    2      NC_007795.1     95777   R       20      32      AAGTCCGTATAATTGTGTAAAAGTAAAAAGGCCATATAACAGTCCTTTTACGGTACAATGTTTTTAACGACAAAAACA
    3      NC_007795.1     132007  L       21      24      TAAAAAAAGCAGGAAGTTTTACCTTCCCACCATAAAATATGAAGAACCTGAAA
    4      NC_007795.1     132011  R       23      27      GGTTCTCCACCAAATGTGGTGGGTATATAATTTAAAGAACAATTTT
    5      NC_007795.1     360608  L       31      34      AAAAGGTGTAAACTCAGTCATAATTATCAATTCCTTTCATAAGAATTAACTGTTAACTAGAGTTTACACCAC
    6      NC_007795.1     360623  R       26      42      GGCAGTGTGTCAATAGAAAAATGAGAATTTCTTAAATTTTTTTAATAAATATTACTAAAAGAGTACAGAAGAT
    7      NC_007795.1     512471  L       25      20      AAAAAGCAGGAAGTTTTACCTTCCCACCATAAAAGATGAAGAACC
    8      NC_007795.1     512472  R       18      25      GTTCTCCACCAAATGTGGTGGGTATATAATTTAAAGAACTATTTT
    9      NC_007795.1     515495  L       22      38      ATNTTTAAACAATCAAAAGTGTACATTATTAAATTATCATTTCCA
    ...

Here is a description of each column:

    flank_id - An arbitrary identifier used internally by *mustache*.
    contig - The name of the contig where the flank was identified.
    pos - The exact base pair position of where the clipped-ends begin, the presumed start of the insertion. 
          This location is 0-based with respect to the reference genome used.
    orient - The orientation of the flank.
             'R' refers to right, which suggests the flank runs from the 5' to the 3' direction with respect to the reference genome.
                    These are reconstructions of the 5' end of the candidate insertion.
             'L' refers to left, which suggests the flank runs from the 3' to the 5' direction with respect to the reference genome.
                    These are reconstructions of the 3' end of the candidate insertion.
    softclip_count - The number of reads that support the existence of this insertion
    runthrough_count - The number of reads that run through the position without being softclipped.
                       This is typically 50% of the number of softclipped reads.
    consensus_seq - The consensus sequence of the flank insertion at the specified site.

This intermediate file is likely not valuable to users on its own, which is why we need to pair insertion flanks, as in the following command.

## `pairflanks`
This command pairs identified insertion flanks with eachother to represent the 5' and 3' ends of inserted sequence.

You can run the command as

    mustache pairflanks FLANKSFILE BAMFILE GENOME

Where `FLANKSFILE` is the output of the `findflanks` command, `BAMFILE` is the binary sequence alignment file used as input
for `findflanks`, and `GENOME` is the reference genome that `BAMFILE` was aligned to.

One important additional parameter is 

    -maxdr, --max_direct_repeat_length
    
which specifies the maximum distance that insertion flanks can be from each other in order to consider pairing them together.
Since insertions often cause direct repeats at the insertion site, the location of the insertion flanks are often seperated by several base pairs.
The default for this parameter is 20. This means that our algorithm will not identify true insertions that create direct repeats that exceed 20 base pairs.

The `pairflanks` command uses a variety of techniques to pair flanks with each other. It first filters insertion flanks by the
`--max_direct_repeat_length` parameter. It then does pairwise comparisons between all nearby, oppositely-oriented insertion flanks
to determine if they share inverted repeats at their termini, a common feature of many prokaryotic insertion sequences.
It then prioritizes all pairs by the following parameters, in order: 
1. The length of the inverted repeat identified, if any. Pairs with longer shared inverted repeats receive higher priority.
2. The difference in the length of the recovered insertion flanks. 
If the two insertion flanks have very similar lengths, they are given higher priority than pairs with disparate flank lengths.
3. The difference in the number of softclipped reads that support the insertion flanks. 
Pairs that have a similar number of reads supporting each flank are given higher priority.

We realize that these are somewhat arbitrary filters, and we plan to improve upon this pairing prioritization in the future.

## `inferseq-assembly`

## `inferseq-reference`

## `inferseq-overlap`

## `inferseq-database`

## `formatbam`

## `recall`

## `pipeline`

## `extendpairs`
