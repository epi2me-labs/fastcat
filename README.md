# fastcat

`fastcat` summarises FASTQ and BAM data and can create offline quality-control reports.
Use `fastcat qc` for new analyses; it supports read filtering, optional barcode demultiplexing,
alignment metrics, and coverage summaries where applicable.

## Installation

### Conda

A conda package is available through our conda channel:

```bash
mamba create -n fastcat -c conda-forge -c bioconda -c nanoporetech fastcat
```

### Build from source

`fastcat` depends on bundled third-party sources in git submodules (`third_party/htslib`, `third_party/zlib-ng`, and `third_party/argp-standalone`). Purrflate is an optional acceleration submodule for multithreaded fastq.gz decompression.

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

`--recursive` initialises the public dependencies and deliberately skips the private optional Purrflate submodule. If you already cloned without submodules, run the following to initialise the public dependencies:

```bash
git submodule update --init third_party/htslib third_party/zlib-ng third_party/argp-standalone
```

Internal builds that have the Purrflate submodule use it automatically. Set `NOPURRFLATE=1` to explicitly select the portable zlib FASTQ reader, while retaining Fastcat's normal threading support.

## CLI overview

The command-line interface is split into multiple sub-commands with distinct tasks:

```bash
$ ./fastcat --help
Usage: ./fastcat <subcommand> [options]
Subcommands:
  qc           Unified FASTQ/BAM quality-control and reporting mode
  fastq        (Deprecated) FASTQ concatenation/statistics/filtering mode
  bam          (Deprecated) BAM alignment statistics mode
  lint         Standalone FASTQ low-complexity filter
  index        BAM chunk index tools (bamindex)
  json-schema  Write the summary.json JSON Schema to standard output
Try './fastcat <subcommand> --help' for subcommand options.
```

### `fastcat qc`

Use `fastcat qc` to summarise FASTQ or BAM data and, optionally, create an offline HTML report.
fastcat detects a common FASTQ or BAM format from named input files, directories, and
file-of-filenames. Specify `--input-fmt` when reading input data from standard input, or to
require a particular format.

```bash
# FASTQ quality-control report
fastcat qc reads.fastq.gz --html-report --output fastcat-qc

# BAM quality-control report
fastcat qc alignments.bam --html-report --output fastcat-qc
```

The primary QC output of `fastcat qc` is a summary.json file which includes all information
garnered from the input files.

FASTQ and unaligned BAM inputs can be filtered with the read-filtering options and are written to
standard output.
Use `--reads-output PREFIX` to save them as `PREFIX.fastq.gz`, or `PREFIX.bam` with `--bam-out`,
beneath the report directory.
When read filtering is enabled, read-level QC metrics describe reads that pass the filters; the summary also reports counts of records rejected by each filter.
For a BAM with reference sequences, fastcat collects alignment statistics only; it does not filter or write reads.
The `--region` and `--bed` options can be used to select genomic regions from indexed BAMs.
With one coordinate-sorted reference BAM, fastcat includes reference coverage summaries.
Use `--coverage` to additionally write position-level coverage BED.gz files.
For BAM input, `--poly-a` enables poly-A metrics for alignments with suitable data.

Useful options include `--html-report`, `--output DIRECTORY`, `--threads THREADS`, `--sample NAME`,
`--demultiplex`, and `--legacy-outputs`. Run `fastcat qc --help` for the complete option reference.
Use `--jsonl-progress PATH` when a wrapper needs machine-readable lifecycle and progress events;
each line in that file is a complete JSON object. Progress events are emitted at most once per
second, independently of the less-frequent human-readable terminal progress log. Every event has
`schema_version`, a per-run monotonically increasing `sequence`, an ISO 8601 UTC `timestamp`, and
an `event` type; remaining fields are specific to that event.

### Legacy commands

`fastcat fastq` and `fastcat bam` are deprecated compatibility commands.
New QC analyses should use `fastcat qc`; the older commands remain available for read concatenation, filtering, and existing
automation.
They retain their historical text-output default, which can be disabled with `--no-legacy-outputs`.
Use `fastcat fastq --help` or `fastcat bam --help` for their option reference.


## Outputs

`fastcat qc` writes its results under `-o`/`--output DIRECTORY` (default `fastcat-qc`; the directory
must not already exist).
Results are first written to `DIRECTORY.tmp` and published as `DIRECTORY` only after a successful run.
Use `--keep-failed-output` to retain that temporary directory for debugging after a failure.
Input and output failures use the conventional `sysexits` status values where applicable; an
interrupted run exits with `128 +` its signal number after removing its temporary output.
With `--demultiplex`, it writes one output directory for each observed
barcode under `<output>/{unclassified,barcodeXXXX/}`.

Each output directory contains a `summary.json` file, a machine-readable summary of the analysis.
The `--html-report` option adds an offline interactive report.
With `--legacy-outputs` enabled, fastcat also writes legacy aggregate text outputs: `per_file.tsv`,
`runids.tsv`, `basecallers.tsv`, histogram files, BAM `flagstats.tsv`, and coverage summary and
distribution tables where applicable. Their aggregate information is included in `summary.json`.
`per_read.tsv`, requested separately with `--per-read`, is a detailed export and is not controlled
by `--legacy-outputs`.
Enabling the `--coverage` option results in detailed coverage BED.gz traces being written, in the case of a single coordinate-sorted reference BAM.

### Read output

For FASTQ and unaligned BAM input, `fastcat qc` writes filtered reads to standard output in FASTQ
format, or unaligned BAM with `--bam-out`.
It does not write a read stream for reference-bearing BAM input.
With `--demultiplex`, reads are written to barcode-specific output directories instead.

### summary.json

The JSON summary contains statistics about the input data, including read and base yield, read length,
mean quality, filtering outcomes, sequencing-over-time metrics, channel metrics, and quality
diagnostics.
Reference-bearing BAM input also adds alignment and coverage summaries.
Most values are directly readable JSON; larger arrays are encoded compactly for efficient storage.

> The following documentation is generated automatically from a JSON schema document describing the file contents.
> The JSON Schema document can be obtained by running `fastcat json-schema` from the CLI.
> Use this file if you need a machine readable specification of the `summary.json` files.

#### Dataframe encoding

For large array-like data the `summary.json` files use a self-describing base64-encoding.
The descriptions contain the data type, the shape of the data and finally the encoded data.

The table below enumerates various definitions seen in the JSON schema document describing `summary.json` files.
Examples of these structures are given folowing the table.

| Schema definition | Type | Description |
| --- | --- | --- |
| `array` | `object` | Self-describing array container. It has `kind: "array"`, `dtype`, `shape`, and `data`; `name` is optional. Numeric payloads are base64-encoded little-endian values, while string payloads are direct JSON arrays. |
| `dataframe` | `object` | Column-oriented table with `kind: "dataframe"` and a two-dimensional `[rows, columns]` shape. Every column is an independently decodable one-dimensional `array`. |
| `adaptiveColumn` | `object` | Named one-dimensional numeric array used for unsigned integer columns. Its dtype is selected per output as `<u1`, `<u2`, `<u4`, or `<u8`. |
| `floatColumn` | `object` | Named one-dimensional numeric array with fixed little-endian float64 dtype, `<f8`. |
| `stringColumn` | `object` | Named one-dimensional string array with `dtype: "string"` and direct JSON string values. |

For example, 1D and 2D arrays are represented as `array` objects:
```
{
  "one_dimensional": {
    "name": "my_1d_array", (optional)
    "kind": "array",
    "dtype": "<u2",
    "shape": [3],
    "data": "AQACAAMA"
  }
  "two_dimensional": {
    "kind": "array",
    "dtype": "<f4",
    "shape": [2, 2],
    "data": "AACAPwAAAEAAAEBAAACAQA=="
  }
}
```

To optimize filesizes, the data type size, as represented by the `dtype` field in the above examples, may differ across files for the same attribute.
Parsers should therefore not assume that a given array or data frame column is always encoded with the same data type, but rather inspect the `dtype` field.

Tabular, `dataframe` objects, data are simply lists of named `array` objects:

```
{
  table: {
    "kind": "dataframe",
    "shape": [2, 2],
    "columns": [
      {
        "name": "chrom",
        "kind": "array",
        "dtype": "string",
        "shape": [2]
        "data": ["chr1", "chr2"]
      },
      {
        "name": "mean_depth",
        "kind": "array",
        "dtype": "<u4",
        "shape": [2]
        "data": "ZAAAAMgAAAA="
      }
    ]
  }
}
```

<!-- BEGIN GENERATED SUMMARY.JSON REFERENCE -->

#### `run_summary` — Run summary

The principal human-readable overview for this output sample or barcode: accepted-read yield and count, filtering outcomes, sequencing-run metadata, and N/L read-length statistics.

| Field | Type | Description |
| --- | --- | --- |
| `sample_id` | `string` or `null` | User-supplied sample identifier, or null when none was supplied. |
| `input_type` | `string` or `null` | Detected input type, for example fastq, unaligned bam, or aligned bam. |
| `barcode_id` | `string` or `null` | Barcode identifier for this output, or null when it is not applicable. |
| `barcode_alias` | `string` or `null` | User-provided barcode alias, or null when no alias is available. |
| `run_ids` | `array` or `null` | Distinct sequencing run identifiers observed in accepted reads, or null when unavailable. |
| `basecaller_models` | `array` or `null` | Distinct basecaller model identifiers observed in accepted reads, or null when unavailable. |
| `modbase_models` | `array` or `null` | Distinct modified-base model identifiers observed in accepted reads, or null when unavailable. |
| `flow_cell_ids` | `array` or `null` | Distinct flow-cell identifiers observed in accepted reads, or null when unavailable. |
| `device_ids` | `array` or `null` | Distinct sequencing-device identifiers observed in accepted BAM reads, or null when unavailable. |
| `run_date` | `string` or `null` | Earliest observed sequencing start date, or null when timestamps are unavailable. |
| `run_duration_s` | `number` | Elapsed sequencing duration in seconds; -1 denotes multiple flow cells and 0 is used when no valid timing range is available. |
| `read_count` | `nonNegativeInt` | Number of reads that passed fastcat filtering. |
| `base_count` | `nonNegativeInt` | Total bases in reads that passed fastcat filtering. |
| `records_seen` | `nonNegativeInt` | Number of input records considered before read-level filtering. |
| `records_passed_filters` | `nonNegativeInt` | Number of records that passed all read-level filters. |
| `records_too_long` | `nonNegativeInt` | Number of records rejected by the maximum read-length filter. |
| `records_too_short` | `nonNegativeInt` | Number of records rejected by the minimum read-length filter. |
| `records_low_quality` | `nonNegativeInt` | Number of records rejected by the mean-quality filter or missing usable quality. |
| `records_dust_masked` | `nonNegativeInt` | Number of records rejected by low-complexity (DUST) filtering. |
| `n5_bp` | `nonNegativeInt` | N5 read length in base pairs: 5% of total accepted bases are in reads of at least this length. |
| `n10_bp` | `nonNegativeInt` | N10 read length in base pairs: 10% of total accepted bases are in reads of at least this length. |
| `n25_bp` | `nonNegativeInt` | N25 read length in base pairs: 25% of total accepted bases are in reads of at least this length. |
| `n50_bp` | `nonNegativeInt` | N50 read length in base pairs: 50% of total accepted bases are in reads of at least this length. |
| `n75_bp` | `nonNegativeInt` | N75 read length in base pairs: 75% of total accepted bases are in reads of at least this length. |
| `n90_bp` | `nonNegativeInt` | N90 read length in base pairs: 90% of total accepted bases are in reads of at least this length. |
| `n95_bp` | `nonNegativeInt` | N95 read length in base pairs: 95% of total accepted bases are in reads of at least this length. |
| `l5` | `nonNegativeInt` | Number of longest reads required to reach 5% of accepted bases. |
| `l10` | `nonNegativeInt` | Number of longest reads required to reach 10% of accepted bases. |
| `l25` | `nonNegativeInt` | Number of longest reads required to reach 25% of accepted bases. |
| `l50` | `nonNegativeInt` | Number of longest reads required to reach 50% of accepted bases. |
| `l75` | `nonNegativeInt` | Number of longest reads required to reach 75% of accepted bases. |
| `l90` | `nonNegativeInt` | Number of longest reads required to reach 90% of accepted bases. |
| `l95` | `nonNegativeInt` | Number of longest reads required to reach 95% of accepted bases. |

#### `sample` — Sample metadata

Sample and barcode identifiers for this output.

| Field | Type | Description |
| --- | --- | --- |
| `sample_id` | `string` or `null` | User-supplied sample identifier, or null when none was supplied. |
| `input_type` | `string` or `null` | Detected input type, or null when unavailable. |
| `barcode_id` | `string` or `null` | Barcode identifier for this output, or null when it is not applicable. |
| `barcode_alias` | `string` or `null` | User-provided barcode alias, or null when unavailable. |

#### `files` — Input file summaries

Per-input-file accepted-read summaries. In a demultiplexed run, each row is scoped to the barcode represented by this summary.json file.

| Field | Type | Description |
| --- | --- | --- |
| `path` | `string` | Input-file path as supplied to fastcat. |
| `accepted_read_count` | `nonNegativeInt` | Number of reads from this file that passed fastcat filtering. |
| `accepted_base_count` | `nonNegativeInt` | Total bases in reads from this file that passed fastcat filtering. |
| `min_length_bp` | `nonNegativeInt` | Shortest accepted read length in base pairs. |
| `max_length_bp` | `nonNegativeInt` | Longest accepted read length in base pairs. |
| `mean_quality` | `number` | Arithmetic mean of accepted reads' mean Phred qualities. |
| `run_id_counts` | `counter` | Exact run-ID counts among accepted reads from this file. |
| `basecaller_counts` | `counter` | Exact basecaller-model counts among accepted reads from this file. |

#### `counts` — Aggregate count summaries

Machine-readable aggregate accepted-read and base totals, with frequency counters that preserve the number of reads associated with each observed run ID and basecaller model.

| Field | Type | Description |
| --- | --- | --- |
| `read_count` | `nonNegativeInt` | Number of reads that passed fastcat filtering. |
| `base_count` | `nonNegativeInt` | Total bases in reads that passed fastcat filtering. |
| `run_ids` | `counter` | Exact accepted-read counts by run identifier. |
| `basecallers` | `counter` | Exact accepted-read counts by basecaller model identifier. |

#### `channels` — Per-channel metrics

Sparse channel-and-mux aggregates used for the channel-occupancy maps: accepted read and base totals, plus mean quality and, for aligned input, mean alignment concordance. This value is null when usable channel metadata is unavailable or when an output contains multiple flowcells.

| Field | Type | Description |
| --- | --- | --- |
| `device` | one of `minion`, `flongle`, `promethion`, `unknown` | Flow-cell device family inferred from observed channel numbers. |
| `skipped` | `channelSkipped` | Counts of reads omitted from channel aggregation because channel or mux metadata was unusable. |
| `data` | `dataframe` | Sparse channel dataframe, indexed by physical channel and mux. |

The **`data`** dataframe columns are:

| Field | Type | Description |
| --- | --- | --- |
| `channel` | `adaptiveColumn` | Physical sequencing channel number. |
| `mux` | `adaptiveColumn` | Mux value associated with the channel observation. |
| `read_count` | `adaptiveColumn` | Number of accepted reads assigned to this channel and mux. |
| `base_count` | `adaptiveColumn` | Total accepted bases assigned to this channel and mux. |
| `quality_read_count` | `adaptiveColumn` | Number of reads with a usable mean quality for this channel and mux. |
| `mean_quality` | `floatColumn` | Arithmetic mean read quality for reads counted by quality_read_count. |
| `concordance_read_count` | `adaptiveColumn` | Number of aligned reads with usable concordance for this channel and mux. |
| `mean_concordance_qscore` | `floatColumn` | Arithmetic mean alignment concordance on the Phred scale for contributing reads. |


#### `distributions` — Read statistic histograms

Histogram data for canonical reads: primary mapped and primary unmapped BAM records share the same read-length and mean-quality distributions.

| Field | Type | Description |
| --- | --- | --- |
| `read_length` | optional value (`null` when unavailable) | Read-length histogram with a reverse-cumulative yield trace, or null when no reads are available. |
| `quality` | optional value (`null` when unavailable) | Mean-read-quality histogram, or null when no reads are available. |

The **`read_length.data`** dataframe columns are:

| Field | Type | Description |
| --- | --- | --- |
| `start` | `adaptiveColumn` | Inclusive lower edge of the read-length bin, in base pairs. |
| `end` | `adaptiveColumn` | Exclusive upper edge of the read-length bin, in base pairs; zero denotes an unbounded final bin. |
| `count` | `adaptiveColumn` | Number of reads in the bin. |
| `yield_at_or_above_bp` | `adaptiveColumn` | Total bases in reads with length at or above this bin's lower edge. |


The **`quality.data`** dataframe columns are:

| Field | Type | Description |
| --- | --- | --- |
| `start` | `floatColumn` | Inclusive lower edge of the histogram bin. |
| `end` | `floatColumn` | Exclusive upper edge of the histogram bin. |
| `count` | `adaptiveColumn` | Number of observations in the bin. |


#### `time` — Sequencing metrics over time

Chronological per-bin read and base yield, cumulative totals, and read-statistic quantiles. These values underpin the sequencing-over-time plots. This value is null when usable timing metadata is unavailable.

| Field | Type | Description |
| --- | --- | --- |
| `skipped` | `timeSkipped` | Counts of records excluded from time aggregation because timing metadata was missing or invalid. |
| `data` | `dataframe` | Time dataframe with one row per elapsed-time bin. |

The **`data`** dataframe columns are:

| Field | Type | Description |
| --- | --- | --- |
| `time_bin_start_s` | `floatColumn` | Inclusive elapsed-time start of the bin, in seconds. |
| `time_bin_end_s` | `floatColumn` | Exclusive elapsed-time end of the bin, in seconds. |
| `read_count` | `adaptiveColumn` | Accepted reads whose sequencing start falls in this bin. |
| `base_count` | `adaptiveColumn` | Bases in accepted reads whose sequencing start falls in this bin. |
| `cumulative_read_count` | `adaptiveColumn` | Accepted reads observed from run start through this bin. |
| `cumulative_base_count` | `adaptiveColumn` | Accepted bases observed from run start through this bin. |
| `length_mean` | `floatColumn` | Arithmetic mean read length for reads in this bin, in base pairs. |
| `length_q05` | `floatColumn` | Q05 read-length quantile for reads in this bin, in base pairs. |
| `length_q10` | `floatColumn` | Q10 read-length quantile for reads in this bin, in base pairs. |
| `length_q25` | `floatColumn` | Q25 read-length quantile for reads in this bin, in base pairs. |
| `length_q50` | `floatColumn` | Q50 read-length quantile for reads in this bin, in base pairs. |
| `length_q75` | `floatColumn` | Q75 read-length quantile for reads in this bin, in base pairs. |
| `length_q90` | `floatColumn` | Q90 read-length quantile for reads in this bin, in base pairs. |
| `length_q95` | `floatColumn` | Q95 read-length quantile for reads in this bin, in base pairs. |
| `quality_mean` | `floatColumn` | Arithmetic mean read quality for reads in this bin, on the Phred scale. |
| `quality_q05` | `floatColumn` | Q05 read-quality quantile for reads in this bin, on the Phred scale. |
| `quality_q10` | `floatColumn` | Q10 read-quality quantile for reads in this bin, on the Phred scale. |
| `quality_q25` | `floatColumn` | Q25 read-quality quantile for reads in this bin, on the Phred scale. |
| `quality_q50` | `floatColumn` | Q50 read-quality quantile for reads in this bin, on the Phred scale. |
| `quality_q75` | `floatColumn` | Q75 read-quality quantile for reads in this bin, on the Phred scale. |
| `quality_q90` | `floatColumn` | Q90 read-quality quantile for reads in this bin, on the Phred scale. |
| `quality_q95` | `floatColumn` | Q95 read-quality quantile for reads in this bin, on the Phred scale. |
| `translocation_speed_mean` | `floatColumn` | Arithmetic mean translocation speed for reads in this bin, in bases per second. |
| `translocation_speed_q05` | `floatColumn` | Q05 translocation-speed quantile for reads in this bin, in bases per second. |
| `translocation_speed_q10` | `floatColumn` | Q10 translocation-speed quantile for reads in this bin, in bases per second. |
| `translocation_speed_q25` | `floatColumn` | Q25 translocation-speed quantile for reads in this bin, in bases per second. |
| `translocation_speed_q50` | `floatColumn` | Q50 translocation-speed quantile for reads in this bin, in bases per second. |
| `translocation_speed_q75` | `floatColumn` | Q75 translocation-speed quantile for reads in this bin, in bases per second. |
| `translocation_speed_q90` | `floatColumn` | Q90 translocation-speed quantile for reads in this bin, in bases per second. |
| `translocation_speed_q95` | `floatColumn` | Q95 translocation-speed quantile for reads in this bin, in bases per second. |
| `accuracy_read_count` | `adaptiveColumn` | Aligned reads with usable alignment concordance in this bin. |
| `accuracy_mean` | `floatColumn` | Arithmetic mean alignment concordance for reads in this bin, on the Phred scale. |
| `accuracy_q05` | `floatColumn` | Q05 alignment-concordance quantile for reads in this bin, on the Phred scale. |
| `accuracy_q10` | `floatColumn` | Q10 alignment-concordance quantile for reads in this bin, on the Phred scale. |
| `accuracy_q25` | `floatColumn` | Q25 alignment-concordance quantile for reads in this bin, on the Phred scale. |
| `accuracy_q50` | `floatColumn` | Q50 alignment-concordance quantile for reads in this bin, on the Phred scale. |
| `accuracy_q75` | `floatColumn` | Q75 alignment-concordance quantile for reads in this bin, on the Phred scale. |
| `accuracy_q90` | `floatColumn` | Q90 alignment-concordance quantile for reads in this bin, on the Phred scale. |
| `accuracy_q95` | `floatColumn` | Q95 alignment-concordance quantile for reads in this bin, on the Phred scale. |


#### `heatmaps` — Quality diagnostic heatmaps

Two-dimensional count matrices relating mean read quality to read length, translocation speed, or alignment concordance. Axis metadata identifies the metric and its numeric range.

| Field | Type | Description |
| --- | --- | --- |
| `read_length_quality` | optional value (`null` when unavailable) | Read-length versus mean-quality count matrix. |
| `translocation_speed_quality` | optional value (`null` when unavailable) | Translocation-speed versus mean-quality count matrix. |
| `alignment_accuracy_quality` | optional value (`null` when unavailable) | Alignment-concordance versus mean-quality count matrix; available for aligned BAM data. |

#### `alignment` — Alignment summary metrics

BAM-only record-level mapping counts, canonical-read mapping totals, per-reference flag counts, and distributions for alignment concordance, aligned query coverage, and poly-A length. This section is omitted for FASTQ summaries.

| Field | Type | Description |
| --- | --- | --- |
| `record_counts` | `alignmentCounts` | Aggregate BAM record counts by flag-derived category. |
| `mapped_record_count` | `nonNegativeInt` | Number of mapped BAM records. |
| `mapped_record_percent` | `percentage` | Percentage of all BAM records that are mapped. |
| `primary_mapped_record_count` | `nonNegativeInt` | Number of records that are both primary and mapped. |
| `primary_mapped_record_percent` | `percentage` | Percentage of all BAM records that are both primary and mapped. |
| `read_count` | `nonNegativeInt` | Number of canonical reads represented by the alignment records. |
| `mapped_read_count` | `nonNegativeInt` | Number of canonical reads with a primary mapped record. |
| `mapped_read_percent` | `percentage` | Percentage of canonical reads with a primary mapped record. |
| `per_reference_counts` | `array` | Flag-derived record counts for every reference sequence. |
| `distributions` | `alignmentDistributions` | Alignment-specific read distributions. |

The **`distributions.accuracy.data` and `distributions.query_coverage.data`** dataframe columns are:

| Field | Type | Description |
| --- | --- | --- |
| `start` | `floatColumn` | Inclusive lower edge of the histogram bin. |
| `end` | `floatColumn` | Exclusive upper edge of the histogram bin. |
| `count` | `adaptiveColumn` | Number of observations in the bin. |


The **`distributions.poly_a_length.data`** dataframe columns are:

| Field | Type | Description |
| --- | --- | --- |
| `start` | `adaptiveColumn` | Inclusive lower edge of the read-length bin, in base pairs. |
| `end` | `adaptiveColumn` | Exclusive upper edge of the read-length bin, in base pairs; zero denotes an unbounded final bin. |
| `count` | `adaptiveColumn` | Number of reads in the bin. |


#### `coverage` — Reference coverage

Reference-depth tables for a single coordinate-sorted BAM with reference sequences. Target sets may cover the whole reference, fixed-length segments, or named BED intervals; the summary and cumulative-distribution dataframes directly mirror fastcat's legacy coverage tables. This section is omitted for FASTQ, unaligned BAM, multiple BAM inputs, and unsorted aligned BAM.

| Field | Type | Description |
| --- | --- | --- |
| `targets` | `array` | Coverage tables for all configured global, segment, and BED target sets. |

The **`targets[].summary.data`** dataframe columns are:

| Field | Type | Description |
| --- | --- | --- |
| `chrom` | `stringColumn` | Reference sequence or target-set label; the final aggregate row uses the target-set name. |
| `start` | `adaptiveColumn` | Zero-based inclusive interval start. |
| `end` | `adaptiveColumn` | Zero-based exclusive interval end. |
| `length` | `adaptiveColumn` | Interval length in base pairs. |
| `bases` | `adaptiveColumn` | Total aligned bases contributing depth within the interval. |
| `mean` | `floatColumn` | Mean reference depth across interval positions. |
| `min` | `adaptiveColumn` | Minimum reference depth across interval positions. |
| `max` | `adaptiveColumn` | Maximum reference depth across interval positions. |
| `^[0-9]+x$` | `floatColumn` | Fraction of interval positions covered at or above the named depth threshold. |


The **`targets[].distribution.data`** dataframe columns are:

| Field | Type | Description |
| --- | --- | --- |
| `chrom` | `stringColumn` | Target-set label; coverage distributions are aggregate-only. |
| `depth` | `adaptiveColumn` | Minimum reference depth. |
| `fraction` | `floatColumn` | Fraction of target positions covered at or above depth. |

<!-- END GENERATED SUMMARY.JSON REFERENCE -->




### Histogram outputs

> These outputs are considered legacy outputs and may be removed in future versions. It is recommended to use the `summary.json` output instead.

The program creates a `histograms/` directory (per barcode, when demultiplexed) containing:

* `length.hist` - read length histogram.
* `quality.hist` - read mean base-quality score histogram.
* `accuracy.hist` **(bam only)** - read alignment accuracy histogram.
* `coverage.hist` **(bam only)** - read alignment coverage histogram (the proportion of each read spanned by the alignment and not clipped).
* `polya.hist` **(bam only, with `--poly-a`)** - poly-A tail length histogram.

The format of the histogram files is a tab-separated file of sparse, ordered intervals `[lower, upper)`:

```text
lower    upper    count
```

The final bin may be unbounded, which is signified by a `0` entry for the upper bin edge.

### Legacy per-file outputs

`per_file.tsv` is a legacy per-input processing diagnostic, rather than part of the supported
aggregate JSON interface. It is a tab-separated file with columns:

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
