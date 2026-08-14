# fastcat quality control

Perform sanity checks on sequencing output files, and produce quality control summary reports.

## Introduction

fastcat provides a quick quality-control report for nanopore sequencing reads. It runs locally on
your computer and produces an interactive report that you can open without an internet connection.
EPI2ME workflows run fastcat as part of their own input checks; this standalone tool lets you view
the same report before uploading data. It can therefore identify damaged or unreadable files
locally, rather than after a large upload has started a workflow that cannot continue.

## Choose your input

Select a FASTQ file or a BAM file. fastcat automatically detects the format of a single file.
When you select a directory, choose the expected format so fastcat searches only for that type of file.

- **FASTQ**: use for basecalled reads before alignment. The report includes yield, read length,
  read quality, sequencing-over-time, channel, and quality-diagnostic summaries when that
  information is present in the read headers.
- **BAM**: use for aligned or unaligned reads stored in BAM. For aligned BAM files, the report also includes
  mapping and alignment summaries, alignment concordance, and reference-depth information where
  available.

You can also select a directory containing input files of the chosen format. All selected files
are summarised together.

## Optional settings

- **Demultiplex barcodes**: creates separate summaries for each detected barcode.
- **Search subfolders**: includes matching files in subfolders of a selected input directory.
- **Check input integrity**: stops the analysis if fastcat finds an
  unreadable input file or invalid read records, helping identify issues before data are uploaded
  or used in later analysis.
- **Save reads** keeps a copy of the reads that pass filtering for FASTQ or unaligned BAM
  input. This option is not available for aligned BAM input.

## Your results

The results include an interactive HTML report. It summarises read yield, count, length, quality,
and sequencing performance over time. Where the input contains the relevant information, it also
shows channel occupancy, barcode summaries, alignment performance, and reference coverage.

For definitions of report metrics and further guidance, see the
[fastcat documentation](https://epi2me.nanoporetech.com/epi2me-docs/tools/fastcat/).
