# fastcat

`fastcat` is a single executable with subcommands for FASTQ and BAM summary/statistics workflows.
It was originally written to santize and concatenate .fastq(.gz) files, hence the name.
Contemporary usage is to provide a single tool for summarising and filtering FASTQ and BAM files, with optional demultiplexing by barcode, and production of QC reports.

## Installation

### Conda

A conda package is available through our conda channel:

```bash
mamba create -n fastcat -c conda-forge -c bioconda -c nanoporetech fastcat
```

### Build from source

`fastcat` depends on bundled third-party sources in git submodules (`third_party/htslib` and `third_party/zlib-ng`). Therefore when cloning be sure to also clone these submodules (or initialize them after clone) before building.

Building requires several prerequisites:

- `git` (with submodule support)
- C toolchain: `gcc` (or `clang`), `make`, `bash`, `patch`
- Autotools: `autoconf`, `automake` (provides `autoreconf`/`autoheader`)
- Development libraries: `zlib`, `bzip2`, `xz/lzma`, `libcurl`, `openssl`
- Cmake is required for the optional depenency of `libdeflate`

For example on Debian or Ubuntu systems run:

```bash
sudo apt update
sudo apt install -y git gcc make bash patch autoconf automake \
  zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev
```

Then to clone and build `fastcat` run:

```bash
git clone --recursive https://github.com/epi2me-labs/fastcat.git
cd fastcat
make fastcat
```

If you already cloned without submodules, run the following to initialise them:

```bash
git submodule update --init --recursive
```

## CLI overview

The command-line interface is split into multiple sub-commands with distinct tasks:

```bash
$ ./fastcat --help
Usage: ./fastcat <subcommand> [options]
Subcommands:
  fastq        FASTQ concatenation/statistics/filtering mode (fastcat)
  bam          BAM alignment statistics/filtering mode
  lint         Standalone FASTQ low-complexity filter (fastlint)
  index        BAM chunk index tools (bamindex)
  json-schema  Write the summary.json JSON Schema to standard output
Try './fastcat <subcommand> --help' for subcommand options.
```

## `fastcat fastq / bam` options

`fastcat fastq` and `fastcat bam` are both built on the same statistics engine and therefore share the majority of
their options (output, read filtering, histogram sketch, time-summary and logging options). These common options
are documented once below; each subcommand's own section then covers only the options unique to it. See the
[outputs section](#fastcat-fastq--bam-outputs) for the files they produce.


```text
 General options:
  -e, --force-error          Exit with non-zero status if any files, or records, contained errors.
  -v, --verbose              Verbose output.

 Output options:
  -o, --output=DIRECTORY     Output root directory (default: fastcat-qc). Must not already exist.
  -d, --demultiplex          Separate barcoded samples using fastq header information; outputs go
                             under `<output>/barcodeXXXX/`.
  -s, --sample=SAMPLE NAME   Sample identifier (defaults to `SAMPLE` if files don't provide
                             sample_id metadata).
  -p, --per-read             Write per-read TSV (large). Disabled by default.
      --html-report          Write a single-sample HTML report using aggregate fastcat histograms.
      --no-legacy-outputs    Disable legacy counter and histogram summary outputs; `summary.json`
                             is still written.

 Read filtering options:
  -a, --min-length=MIN READ LENGTH
                             Minimum read length to output (excluded reads remain listed in
                             summaries).
      --max-length=MAX READ LENGTH
                             Maximum read length to output (excluded reads remain listed in
                             summaries). Available as `-b` on `fastq`; `bam` only accepts the long
                             form since `-b` there is `--bed`.
  -q, --min-qscore=MIN READ QSCORE
                             Minimum read Qscore to output (excluded reads remain listed in
                             summaries).
      --dust                 Enable DUST filtering of reads (default: disabled).
      --max-dust=MAX DUST    Maximum proportion of low-complexity regions to allow in reads
                             (default: 0.95).
      --dust-w=DUST W        Window size for DUST filtering (default: 64).
      --dust-t=DUST T        Threshold for DUST filtering (default: 20).
  -g, --read-group=RG        Only process reads from given read group.
      --tag-name=TN          Only process reads with a given tag (see --tag-value).
      --tag-value=VAL        Only process reads with a given tag value.
      --haplotype=VAL        Only process reads from a given haplotype. Equivalent to
                             --tag-name HP --tag-value VAL.
      --recalc-qual          Force recomputing mean quality, else use the 'qs' tag if present.

 Histogram sketch options (control `summary.json` heatmap detail vs smoothness):
      --length-qscore-k=K    Sketch size for length modeling per qscore bin (higher => more
                             memory, less sketch noise; default: 400).
      --length-qscore-qgrid=N
                             Quantile samples converted into each density curve (higher => finer
                             detail, more CPU; default: 257).
      --length-qscore-qsmooth=W
                             Odd moving-average window over quantile slopes (higher =>
                             smoother/less spiky, but flatter peaks; default: 5).

 Time summary options (control `summary.json` time-series resolution and quantile accuracy):
      --time-bin-size=SECONDS
                             Width of each time bin in seconds (higher => coarser time plots,
                             lower => more rows; default: 3600).
      --time-bin-k=K         Sketch size for per-bin length/qscore quantiles (higher => more
                             memory, less sketch noise; default: 200).

 Logging options:
      --log-file=PATH        Write logs to PATH instead of stderr.
      --log-level=LEVEL      Log level: error, warn, info, debug (default: info).
      --log-mono             Force monochrome logging output (disable ANSI colors).

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version
```

### Additional fastcat fastq options

In addition to the [common options](#fastcat-fastq--bam-options) above, `fastq` accepts:

```text
 General options:
      --fofn                 Treat each input as a file-of-filenames source (use '-' to read that
                             list from stdin).
      --input-fmt=fastq|bam  Input format for fastq mode (default: fastq). Use 'bam' for unaligned
                             BAM input.
  -t, --threads=THREADS      Number of threads for output compression (only with --bam-out).
  -x, --recurse              Search directories recursively for files matching the selected input
                             format.

 Output options:
  -B, --bam-out              Output data as unaligned BAM.
  -c, --reads-per-file=NUM   Split reads into files with a set number of reads (default: single
                             file).
```

Reads input data from stdin when `-` is provided as an input path. Also accepts directories as
input and looks for files matching the selected input format in the top-level directory, recursing
into sub-directories when `-x` is given. The command exits non-zero if any file encountered cannot
be read.

### Additional fastcat bam options

In addition to the [common options](#fastcat-fastq--bam-options) above, `bam` accepts:

```text
 General options:
  -r, --region=chr:start-end Genomic region to process.
  -b, --bed=BEDFILE          BED file for regions to process (mutually exclusive with --region).
  -t, --threads=THREADS      Number of threads (for BAM decompression and BED output compression).

 Coverage calculation options:
      --coverage             Enable coverage calculations.
      --coverage-beds=BEDFILE ...
                             BED file(s) for calculating coverage (space-separated list).
      --coverage-names=NAME ...
                             Name(s) for coverage BED file(s) (space-separated list).
      --thresholds=VALUES ...
                             Coverage thresholds to produce sparse cumulative distribution
                             (space-separated integers).
      --segments=LENGTH ...  Segment length(s) for which to produce outputs (space-separated
                             integers).

 Poly-A options:
      --poly-a               Enable poly-A tail length histogram.
      --poly-a-cover=PCT_COVERAGE
                             Reference alignment coverage for acceptance of read (default: 95).
      --poly-a-qual=QUAL     Read mean Q score for acceptance of read (default: 10).
      --poly-a-rev           Allow reverse alignments (useful for cDNA; default is direct RNA).
```


### `fastcat fastq / bam` outputs

Both `fastcat fastq` and `fastcat bam` write their outputs under `-o`/`--output DIRECTORY` (default `fastcat-qc`; must not already exist). If `--demultiplex` is given, the outputs described below are produced once per observed barcode under `<output>/{unclassified,barcodeXXXX/}` instead of directly under `<output>/`.

Every output directory always contains a `summary.json` file (see below).
Unless `--no-legacy-outputs` is passed, the additional legacy TSV/histogram outputs described in this section are also written; pass `--html-report` to additionally write a single-sample HTML QC report.
Outputs are common to both subcommands except where marked **(bam only)**.

### Read output (stdout) — fastq only

By default `fastcat fastq` writes the (filtered) input sequences to `stdout` in FASTQ format (or as unaligned BAM via `-B`/`--bam-out`), which can be recompressed with `gzip` (or more usefully `bgzip`).
If `--demultiplex` is given, output to `stdout` is not produced; reads are instead written under the per-barcode directories.

### summary.json

Every output directory contains a `summary.json` structured summary file.
This contains various statistics about the input data, including the number of reads, total bases, mean read length, mean quality score, and counts of reads filtered by length or quality thresholds.
It contains also all the data from the legacy TSV outputs described below.
The JSON Schema for `summary.json` can be obtained at any time via `./fastcat json-schema`; which contains descriptions of all entries of the file.


### Histogram outputs

> These outputs are considered legacy outputs and may be removed in future versions. It is recommended to use the `summary.json` output instead.

The program creates a `histograms/` directory (per barcode, when demultiplexed) containing:

* `length.hist` - read length histogram.
* `quality.hist` - read mean base-quality score histogram.
* `accuracy.hist` **(bam only)** - read alignment accuracy histogram.
* `coverage.hist` **(bam only)** - read alignment coverage histogram (the proportion of each read spanned by the alignment and not clipped).
* `length.unmap.hist` / `quality.unmap.hist` **(bam only)** - as `length.hist`/`quality.hist` but for the subset of reads that are unmapped; written alongside the mapped-read histograms whenever the input BAM contains unmapped reads.
* `polya.hist` **(bam only, with `--poly-a`)** - poly-A tail length histogram.

The format of the histogram files is a tab-separated file of sparse, ordered intervals `[lower, upper)`:

```text
lower    upper    count
```

The final bin may be unbounded, which is signified by a `0` entry for the upper bin edge.

### Per-file outputs

The `per_file.tsv` is a tab-separated file with columns:

```text
filename        n_seqs  n_bases min_length      max_length      mean_quality    f_file_ok       f_stream_error  f_qual_missing  f_qual_truncated        f_unknown_error r_record_seen   r_record_ok     r_too_long r_too_short     r_low_quality   r_dust_masked
test/data/bc0.fastq.gz  10      4599    191     965     9.86    1       0       0       0       0       10      10      0       0       0       0
...
```
where the `mean_quality` column is the mean of the per-read `mean_quality` values.
The `f_*` fields detail file-level parse errors.
`r_record_seen` is records parsed/considered, while `r_record_ok` is records that passed filters and were used downstream.

The `basecallers.tsv` and `runids.tsv` files contain per-file counts of unique basecallers and MinKNOW run identifiers observed in the input files.
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

The (optional, enabled with `-p`/`--per-read`) `per_read.tsv` is a tab-separated file with columns:

```text
read_id filename        runid   read_length     mean_quality    channel read_number     start_time
32e13a1c-4171-4706-b6ce-a32c0f65fa16    test/data/bc0.fastq.gz  5a21d8a6996146deceeaea3784244c52741cae93        326     12.31   282     9       2021-04-20T17:00:40Z
b87f011e-b802-4993-8f56-fd240b2e784f    test/data/bc0.fastq.gz  5a21d8a6996146deceeaea3784244c52741cae93        407     9.13    213     19      2021-04-20T17:00:41Z
6f64aedb-bb8e-4777-b494-43e661841e06    test/data/bc0.fastq.gz  5a21d8a6996146deceeaea3784244c52741cae93        355     9.98    67      13      2021-04-20T17:00:41Z
c372fb2c-dd45-4feb-81b2-c167c3d1ce93    test/data/bc0.fastq.gz  5a21d8a6996146deceeaea3784244c52741cae93        317     8.14    337     18      2021-04-20T17:00:41Z
18d04e8d-2816-4986-8e1b-e5be676837fc    test/data/bc0.fastq.gz  5a21d8a6996146deceeaea3784244c52741cae93        965     7.57    507     18      2021-04-20T17:00:41Z
...
```

**(bam only)** additional alignment columns are appended: `ref`, `coverage`, `ref_coverage`, `qstart`, `qend`, `rstart`, `rend`, `aligned_ref_len`, `direction`, `length`, `match`, `ins`, `del`, `sub`, `iden`, `acc`, `duplex`.

The mean quality is defined as:
```
-10 * log10(mean(10^(Q/-10)))
```

where `Q` are the set of all per-base quality scores for the read.
In many cases (where relevant metadata is found in the read header) the mean quality will not be recomputed but taken from the basecaller's value, which may be
subtly different.
Recalculation of the mean quality can be forced with `--recalc-qual`.

### Alignment flag outputs — bam only

> This output is considered a legacy output and may be removed in future versions. It is recommended to use the `summary.json` output instead.

The file `flagstats.tsv` provides similar output to `samtools flagstats`, enumerating counts for different SAM flags for each reference sequence:

```text
ref     sample_name     total   primary secondary       supplementary   unmapped        qcfail  duplicate       duplex  duplex_forming
cneofor_1       0       0       0       0       0       0       0       0       0
cneofor_2       0       0       0       0       0       0       0       0       0
cneofor_3       1       0       0       1       0       0       0       0       0
ecoli1          373     347     5       21      0       0       0       0       0
ecoli2          0       0       0       0       0       0       0       0       0
```

### Alignment coverage outputs — bam only

> The JSON summary contains the contents of the `{name}.summary.txt` and `{name}.dist.txt` files referenced below.
> It does not contain the complete contents of the BED formatted files.

The program can optionally output coverage information (`--coverage`, provided the input BAM is coordinate sorted).
The outputs are inspired by those from [mosdepth](https://github.com/brentp/mosdepth) but have been adapted to be more readily interpretable in line with user expectations.
A set of files is produced for each user-provided BED file (`--coverage-beds`/`--coverage-names`), as well as for the whole genome. The files are:

* `{name}.bed.gz` - a bgzip-compressed BED file with position-level accuracy.
* `{name}.fwd.bed.gz` - as above but only for forward strand alignments.
* `{name}.rev.bed.gz` - as above but only for reverse strand alignments.
* `{name}.summary.txt` - a tab-separated summary of coverage for each region in
  the BED file. And a final total line for all regions in the BED.
* `{name}.dist.txt` - a tab-separated file listing the cumulative distribution function
  of coverage across all regions in the BED file (i.e. "what proportion of positions
  are covered at at least X depth?").

The summary file contains user-configurable entries (`--thresholds`) enumerating a sparse coverage CDF:
the default is to output the proportion of positions covered at 1x, 5x, 10x, 20x, 30x, and 40x.
This is done independently for all regions in the BED, as well as the final total line.


## fastcat lint

`fastcat lint` is a lightweight FASTQ filter for removing reads dominated by low-complexity sequence. It applies
the SDUST algorithm to each read, writes reads at or below the configured masked-base proportion to standard output,
and reports discarded reads on standard error. It is useful as a small pipeline stage when low-complexity reads need
to be excluded before downstream processing.

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

## fastcat index

`fastcat index` splits an unaligned, unsorted BAM logically into record-count chunks without rewriting the source
file. `build` creates a `.bci` index of BGZF offsets, `fetch` streams one indexed chunk as BAM to standard output,
and `dump` makes the index contents inspectable as text. This is intended for scattering record-based BAM work
across parallel workflow tasks; it is not a coordinate index and does not keep records for a query together.

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
