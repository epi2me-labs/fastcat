# fastcat

A set of simply utilities for creating summaries from standard bioinformatics formats.

### Installation

All tools are distributed in a single package from our conda channel, they can be installed
into an isolated conda environment with:

```
mamba create -n fastcat -c conda-forge -c bioconda -c nanoporetech fastcat
```

#### Compilation

Although not recommended, compilation from source is via make:

```
make fastcat bamstats
```

Several libraries are assumed to be present on the system for linking.

### fastcat

This eponymous tool concatenates .fastq(.gz) files whilst creating a summary
of the sequences. Can also demultiplex reads according to Guppy/MinKNOW
.fastq record headers.

```
Usage: fastcat [OPTION...]
            reads1.fastq(.gz) reads2.fastq(.gz) dir-with-fastq ...
fastcat -- concatenate and summarise .fastq(.gz) files.

  -a, --min_length=MIN READ LENGTH
                             minimum read length to output (excluded reads
                             remain listed in summaries).
  -b, --max_length=MAX READ LENGTH
                             maximum read length to output (excluded reads
                             remain listed in summaries).
  -d, --demultiplex=OUT DIR  Separate barcoded samples using fastq header
                             information. Option value is top-level output
                             directory.
  -f, --file=FILE SUMMARY    Per-file summary output
      --histograms=DIRECTORY Directory for outputting histogram information.
                             When --demultiplex is enabled histograms are
                             written to per-sample demultiplexed output
                             directories. (default: fastcat-histograms)
  -H, --reheader             Rewrite fastq header comments as SAM tags (useful
                             for passing through minimap2).
  -q, --min_qscore=MIN READ QSCOROE
                             minimum read Qscore to output (excluded reads
                             remain listed in summaries).
  -r, --read=READ SUMMARY    Per-read summary output
  -s, --sample=SAMPLE NAME   Sample name (if given, adds a 'sample_name'
                             column).
  -x, --recurse              Search directories recursively for '.fastq',
                             '.fq', '.fastq.gz', and '.fq.gz' files.
  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version

Mandatory or optional arguments to long options are also mandatory or optional
for any corresponding short options.

Input files may be given on stdin by specifing the input as '-'. Also accepts
directories as input and looks for .fastq(.gz) files in the top-level
directory. Recurses into sub-directories when the -x option is given. The
command will exit non-zero if any file encountered cannot be read.
```

The program writes the input sequences to `stdout` in .fastq format to be
recompressed with `gzip` (or more usefully `bgzip`).

The `per-read.txt` is a tab-separated file with columns:

```
read_id        filename                read_length  mean_quality
SRR12447496.1  SRR12447496_1.fastq.gz  531          14.03
SRR12447496.2  SRR12447496_1.fastq.gz  513          13.91
SRR12447496.3  SRR12447496_1.fastq.gz  473          14.70
...
```

The mean quality is defined as:
```
-10 * log10(mean(10^(Q/-10)))
```

where `Q` are the set of all per-base quality scores for the read.

The `per-file.txt` is also a tab-separated file with columns:

```
filename                n_seqs  n_bases  min_length  max_length  mean_quality
SRR12447496_1.fastq.gz  16048   8090160  434         697         13.10
SRR12447498_1.fastq.gz  16203   8049713  421         697         13.25
SRR12447499_1.fastq.gz  15484   7812439  424         612         13.16
...
```
where the `mean_quality` column is the mean of the per-read `mean_quality` values.

Additionally as its a common thing to want to do, the program will write
the two files:

* `length.hist` - read length histogram, and
* `quality.hist` - read mean base-quality score histogram.

When data is demultiplexed one such file will be written to the demultiplexed
samples' directories. When demultiplexing is not enabled the files will be
placed in a directory according to the `--histograms` option. The format of the
histogram files is a tab-separated file of sparse, ordered intervals `[upper, lower)`:

```
lower    upper    count
```

The final bin may be unbounded, which is signified by a `0` entry for the upper
bin edge.


### bamstats

The `bamstats` utility is a re-implementation of the `stats_from_bam` program
from [pomoxis](github.com/nanoporetech/pomoxis). It creates read-level summary
statistics of alignments found in a BAM file and reports these in a TSV file.

Additionally as its a common thing to want to do, the program will write
the four files:

* `length.hist` - read length histogram,
* `quality.hist` - read mean base-quality score histogram,
* `accuracy.hist` - read alignment accuracy histogram, and
* `coverage.hist` - read alignment coverage histogram.

These files are as described for the `fastcat` program.

```
Usage: bamstats [OPTION...] <reads.bam>
bamstats -- summarise rears/alignments in one or more BAM files.

 General options:
  -f, --flagstats=FLAGSTATS  File for outputting alignment flag counts.
      --histograms=DIRECTORY Directory for outputting histogram information.
                             (default: bamstats-histograms)
  -r, --region=chr:start-end Genomic region to process.
  -s, --sample=SAMPLE NAME   Sample name (if given, adds a 'sample_name'
                             column).
  -t, --threads=THREADS      Number of threads for BAM processing.

 Read filtering options:

  -g, --read_group=RG        Only process reads from given read group.
      --haplotype=VAL        Only process reads from a given haplotype.
                             Equivalent to --tag_name HP --tag_value VAL.
      --tag_name=TN          Only process reads with a given tag (see
                             --tag_value).
      --tag_value=VAL        Only process reads with a given tag value.
  -u, --unmapped             Include unmapped/unplaced reads in output.

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version

Mandatory or optional arguments to long options are also mandatory or optional
for any corresponding short options.

The program creates a simple TSV file containing statistics for each primary
alignment stored within the input BAM files.
```


### bamindex

The `bamindex` program is a rather curious program that will create a positional index
of alignments in a BAM file. It is intended to be used  within workflows to allow
parallel processing of records in a BAM file, each worker processing a contiguous chunk
of the file. This is most useful with unaligned BAM files.

The program was insired by [bri](https://github.com/jts/bri) by Jared Simpson at [OICR](https://oicr.on.ca/);
which is far cooler.

There are three subcommands:

**bamindex index**

```
$ ./bamindex build --help
Usage: build [OPTION...] <reads.bam>
bamindex build -- create a BAM index corresponding to batches of records.

 General options:
  -c, --chunk_size=SIZE      Number of records in a chunk.
  -t, --threads=THREADS      Number of threads for BAM processing.

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version

Mandatory or optional arguments to long options are also mandatory or optional
for any corresponding short options.

The program creates a simple index of file offsets for each of every (n * M)th
alignment record. No care is taken to keep records corresponding to the same
query together, or any other such niceities. Its intended to be used simply
with unaligned, unsorted BAMs.
```

**bamindex fetch**

```
$ ./bamindex fetch --help
Usage: fetch [OPTION...] <reads.bam.bci>
bamindex fetch -- fetch records from a BAM according to an index.

 General options:
  -c, --chunk=SIZE           Chunk index to retrieve.
  -t, --threads=THREADS      Number of threads for BAM processing.

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version

Mandatory or optional arguments to long options are also mandatory or optional
for any corresponding short options.

The program simply will fetch a batch of records from a BAM fileusing and index
and a chunk ID.
```

**bamindex dump**

```
$ ./bamindex dump --help 
Usage: dump [OPTION...] <reads.bam.bci>
bamindex dump -- dump a BAM chunk index to stdout as text.

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version

The program simply writes the contents of an index to stdout for human
inspection. It has no other purpose.
```
