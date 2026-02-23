# fastcat

`fastcat` is a single executable with subcommands for FASTQ and BAM summary/statistics workflows.
It was originally written to santize and concatenate .fastq(.gz) files, hence the name.

## Installation

### Conda

```bash
mamba create -n fastcat -c conda-forge -c bioconda -c nanoporetech fastcat
```

### Build from source
When building from source with:

```bash
make fastcat
```

several libraries are expected to be present on the system.
Enquire within the `Makefile` for more details.

## CLI overview

The command-line interface is split into multiple sub-commands with distinct tasks:

```bash
$ ./fastcat --help
Usage: ./fastcat <subcommand> [options]
Subcommands:
  fastq   FASTQ concatenation/statistics/filtering mode (fastcat)
  bam     BAM stats/filtering mode (bamstats)
  lint    Standalone FASTQ low-complexity filter (fastlint)
  index   BAM chunk index tools (bamindex)
  cover   BAM reference coverage summaries (bamcoverage)
Try './fastcat <subcommand> --help' for subcommand options.
```

## fastcat fastq

Concatenate and optionally filter FASTQ, optionally demultiplex by barcode, and write summary outputs.

```bash
$ ./fastcat fastq --help
Usage: fastq [OPTION...]
            reads1.fastq(.gz) reads2.fastq(.gz) dir-with-fastq ...
fastcat -- concatenate and summarise .fastq(.gz) files.

 Output options:

 General options:
  -B, --bam-out              Output data as unaligned BAM.
  -c, --reads-per-file=NUM   Split reads into files with a set number of reads
                             (default: single file).
  -d, --demultiplex          Separate barcoded samples using fastq header
                             information; outputs go under
                             <output>/barcodeXXXX/.
  -e, --force-error          Exit with non-zero status if any files, or
                             records, contained errors.
      --fofn                 Treat each input as a file-of-filenames source
                             (use '-' to read that list from stdin).
  -o, --output=DIRECTORY     Output root directory (default: fastcat-qc). Must
                             not already exist.
  -p, --per-read             Write per-read TSV (large). Disabled by default.
  -s, --sample=SAMPLE NAME   Sample name (if given, adds a 'sample_name'
                             column).
  -t, --threads=THREADS      Number of threads for output compression (only
                             with --bam-out).
  -v, --verbose              Verbose output.
  -x, --recurse              Search directories recursively for '.fastq',
                             '.fq', '.fastq.gz', and '.fq.gz' files.

 Read filtering options:
  -a, --min-length=MIN READ LENGTH
                             Minimum read length to output (excluded reads
                             remain listed in summaries).
  -b, --max-length=MAX READ LENGTH
                             Maximum read length to output (excluded reads
                             remain listed in summaries).
      --dust                 Enable DUST filtering of reads (default:
                             disabled).
      --dust-t=DUST T        Threshold for DUST filtering (default: 20).
      --dust-w=DUST W        Window size for DUST filtering (default: 64).
  -g, --read-group=RG        Only process reads from given read group.
      --haplotype=VAL        Only process reads from a given haplotype.
                             Equivalent to --tag-name HP --tag-value VAL.
      --max-dust=MAX DUST    Maximum proportion of low-complexity regions to
                             allow in reads (default: 0.95).
  -q, --min-qscore=MIN READ QSCORE
                             Minimum read Qscore to output (excluded reads
                             remain listed in summaries).
      --recalc-qual          Force recomputing mean quality, else use 'qs' tag
                             if present.
      --tag-name=TN          Only process reads with a given tag (see
                             --tag-value).
      --tag-value=VAL        Only process reads with a given tag value.

 Logging options:
      --log-file=PATH        Write logs to PATH instead of stderr.
      --log-level=LEVEL      Log level: error, warn, info, debug (default:
                             info).
      --log-mono             Force monochrome logging output (disable ANSI
                             colors).

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version
```


The program writes the input sequences to `stdout` by default in .fastq format, which can be recompressed with `gzip`
(or more usefully `bgzip`).
Ancilliary summary outputs are written to the `<output>` directory; these are described below.
If the `--demultiplex` option is provided the output to `stdout` is not produced and instead outputs for each observed barcode (sample), is produced in `<output>/{unclassified,barcodeXXXX/}`.
In this case the ancilliary outputs are produced per barcode.

### Histogram outputs

Additionally as its a common thing to want to do, the program will create two files:

* `histograms/length.hist` - read length histogram, and
* `histograms/quality.hist` - read mean base-quality score histogram.

When data is demultiplexed one such file will be written to the demultiplexed
samples' directories. When demultiplexing is not enabled the files will be
placed in a directory according to the `--histograms` option. The format of the
histogram files is a tab-separated file of sparse, ordered intervals `[lower, uppper)`:

```text
lower    upper    count
```

The final bin may be unbounded, which is signified by a `0` entry for the upper
bin edge.

### Per-file outputs

The `per-file.txt` is also a tab-separated file with columns:

```text
filename        n_seqs  n_bases min_length      max_length      mean_quality    f_file_ok       f_stream_error  f_qual_missing  f_qual_truncated        f_unknown_error r_record_ok     r_too_long r_too_short     r_low_quality   r_dust_masked
test/data/bc0.fastq.gz  10      4599    191     965     9.86    1       0       0       0       0       10      0       0       0       0
...
```
where the `mean_quality` column is the mean of the per-read `mean_quality` values.
The `f_*` fields detail various error and filtering read counts.

The `basecallers.txt` and `runids.tsv` files contain per-file counts of unique basecallers and MinKNOW run identifiers observed in the input files.
Their structures are, respectively:

```text
filename        basecaller      count
test/data/bc0.fastq.gz          10
...
```

and

```text
filename        run_id  count
test/data/bc0.fastq.gz  5a21d8a6996146deceeaea3784244c52741cae93        10
```

### Per-read output

The (optional) `per-read.txt` is a tab-separated file with columns:

```text
read_id filename        runid   read_length     mean_quality    channel read_number     start_time
32e13a1c-4171-4706-b6ce-a32c0f65fa16    test/data/bc0.fastq.gz  5a21d8a6996146deceeaea3784244c52741cae93        326     12.31   282     9       2021-04-20T17:00:40Z
b87f011e-b802-4993-8f56-fd240b2e784f    test/data/bc0.fastq.gz  5a21d8a6996146deceeaea3784244c52741cae93        407     9.13    213     19      2021-04-20T17:00:41Z
6f64aedb-bb8e-4777-b494-43e661841e06    test/data/bc0.fastq.gz  5a21d8a6996146deceeaea3784244c52741cae93        355     9.98    67      13      2021-04-20T17:00:41Z
c372fb2c-dd45-4feb-81b2-c167c3d1ce93    test/data/bc0.fastq.gz  5a21d8a6996146deceeaea3784244c52741cae93        317     8.14    337     18      2021-04-20T17:00:41Z
18d04e8d-2816-4986-8e1b-e5be676837fc    test/data/bc0.fastq.gz  5a21d8a6996146deceeaea3784244c52741cae93        965     7.57    507     18      2021-04-20T17:00:41Z
...
```

The mean quality is defined as:
```
-10 * log10(mean(10^(Q/-10)))
```

where `Q` are the set of all per-base quality scores for the read.
In many cases (where relevant metadata is found in the read header) the mean quality will not be recomputed but taken from the basecaller's value, which may be
subtly different.
Recalculation of the mean quality can be forced with the `--recalc-qual`.

## fastcat bam

Compute BAM-derived summaries, filtering, optional coverage, and optional per-read tables.


```bash
$ ./fastcat bam --help
Usage: bam [OPTION...] <reads.bam>
bamstats -- summarise reads/alignments in an input BAM file.

 Output options:

 General options:
  -b, --bed=BEDFILE          BED file for regions to process.
  -d, --demultiplex          Separate outputs per barcode under
                             <output>/barcodeXXXX/.
  -e, --force-error          Exit with non-zero status if parsing/filtering
                             errors occur.
  -o, --output=DIRECTORY     Output root directory (default: fastcat-qc). Must
                             not already exist.
  -p, --per-read             Write per-read TSV (large). Disabled by default.
  -r, --region=chr:start-end Genomic region to process.
  -s, --sample=SAMPLE NAME   Sample name (if given, adds a 'sample_name'
                             column).
  -t, --threads=THREADS      Number of threads (for BAM decompression and BED
                             output compression).
  -v, --verbose              Verbose output.

 Read filtering options:
  -a, --min-length=MIN READ LENGTH
                             Minimum read length to output (excluded reads
                             remain listed in summaries).
      --dust                 Enable DUST filtering of reads (default:
                             disabled).
      --dust-t=DUST T        Threshold for DUST filtering (default: 20).
      --dust-w=DUST W        Window size for DUST filtering (default: 64).
  -g, --read-group=RG        Only process reads from given read group.
      --haplotype=VAL        Only process reads from a given haplotype.
                             Equivalent to --tag-name HP --tag-value VAL.
      --max-dust=MAX DUST    Maximum proportion of low-complexity regions to
                             allow in reads (default: 0.95).
      --max-length=MAX READ LENGTH
                             Maximum read length to output (excluded reads
                             remain listed in summaries).
  -q, --min-qscore=MIN READ QSCORE
                             Minimum read Qscore to output (excluded reads
                             remain listed in summaries).
      --recalc-qual          Force recomputing mean quality, else use 'qs' tag
                             in BAM if present.
      --tag-name=TN          Only process reads with a given tag (see
                             --tag-value).
      --tag-value=VAL        Only process reads with a given tag value.

 Coverage calculation options:
      --coverage             Enable coverage calculations.
      --coverage-beds=BEDFILE ...
                             BED file(s) for calculating coverage
                             (space-separated list).
      --coverage-names=NAME ...   Name(s) for coverage BED file(s)
                             (space-separated list).
      --segments=LENGTH ...  Segment length(s) for which to produce outputs
                             (space-separated integers).
      --thresholds=VALUES ...   Coverage thresholds to produce sparse
                             cumulative distribution (space-separated
                             integers).

 Poly-A options:
      --poly-a               Enable poly-A tail length histogram.
      --poly-a-cover=PCT_COVERAGE
                             Reference alignment coverage for acceptance of
                             read (default: 95).
      --poly-a-qual=QUAL     Read mean Q score for acceptance of read (default:
                             10).
      --poly-a-rev           Allow reverse alignments (useful for cDNA; default
                             is direct RNA).

 Logging options:
      --log-file=PATH        Write logs to PATH instead of stderr.
      --log-level=LEVEL      Log level: error, warn, info, debug (default:
                             info).
      --log-mono             Force monochrome logging output (disable ANSI
                             colors).

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version
```

---

## fastcat lint

Standalone SDUST-based low-complexity FASTQ filtering.

```text
$ ./fastcat lint --help
Usage: lint [OPTION...] <reads.fastq>
fastlint -- apply sdust algorithm to input files.

 General options:
  -p, --max-proportion=PROPORTION
                             Maximum allowable proportion of masked bases in a
                             read to keep the read (default: 0.95).
  -t, --threshold=THRESHOLD  Threshold for repetition (default: 20).
  -w, --window=WINDOW        Window size (default: 64).

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version
```

---

## fastcat cover

Standalone BAM coverage summaries.

```text
$ ./fastcat cover --help
Usage: cover [OPTION...] <reads.bam>
bamcoverage -- summarise coverage statistics from a BAM file.

 General options:
  -b, --beds=BEDFILE ...     BED files for regions (space-separated list).
  -n, --names=NAME...        Names associated with the BED files
                             (space-separated list).
  -t, --threads=THREADS      Number of threads for BAM processing.

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version
```

---

## fastcat index

Chunk index tooling for BAM.

```text
$ ./fastcat index --help
* bamindex help          Print general help or help about a subcommand.
* bamindex build         Build a BAM index.
* bamindex fetch         Fetch records from a BAM using an index.
* bamindex dump          Dump an index fetch to text.
```

### fastcat index build

```text
$ ./fastcat index build --help
Usage: build [OPTION...] <reads.bam>
bamindex build -- create a BAM index corresponding to batches of records.

 General options:
  -c, --chunk-size=SIZE      Number of records in a chunk.
  -t, --threads=THREADS      Number of threads for BAM processing.

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version
```

### fastcat index fetch

```text
$ ./fastcat index fetch --help
Usage: fetch [OPTION...] <reads.bam.bci>
bamindex fetch -- fetch records from a BAM according to an index.

 General options:
  -c, --chunk=SIZE           Chunk index to retrieve.
  -t, --threads=THREADS      Number of threads for BAM processing.

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version
```

### fastcat index dump

```text
$ ./fastcat index dump --help
Usage: dump [OPTION...] <reads.bam.bci>
bamindex dump -- dump a BAM chunk index to stdout as text.

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version
```
