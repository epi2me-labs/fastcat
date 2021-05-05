# fastcat

A simple utility to concatenate .fastq(.gz) files whilst creating a summary
of the sequences:

```
Usage: fastcat [OPTION...] reads1.fastq(.gz) reads2.fastq(.gz) ...
fastcat -- concatenate and summarise .fastq(.gz) files.

  -f, --file=FILE SUMMARY    Per-file summary output
  -r, --read=READ SUMMARY    Per-read summary output
  -s, --sample=SAMPLE NAME   Sample name (if given adds a 'sample_name'
                             column)
  -x, --recurse              Search directories recursively for '.fastq',
                             '.fq', '.fastq.gz', and '.fq.gz' files.
  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version

Mandatory or optional arguments to long options are also mandatory or optional
for any corresponding short options.

Input files may be given on stdin by specifing the input as '-'.  When the -x
option is given inputs may be directories.
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

