# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.18.1]
### Fixed
- 'run_id' instead of 'basecaller' as column name in basecaller summary output header line.
- `(null)` in FASTQ header comments when run with `-H` on files that had `basecall_model_version_id=...` as only header comment.

## [v0.18.0]
### Added
- Basecaller summary information similar to runid summary.
- RNA poly-A tail length histogram output.
### Fixed
- Random output for runid when not found in header.

## [v0.17.1]
### Added
- `--runids` option to `bamstats` for enumerating detected run identifiers.

## [v0.17.0]
### Added
- `--reads_per_file` option can split inputs into batched files when demultiplexing. Users should use Unix `split` with piped output.
- `--runids` option to output a file enumerating detected run identifiers.
### Changed
- Per-file read statistics now relate to filtered reads only.
- Link `fastcat` against zlib-ng for an even faster cat.

## [v0.16.8]
### Changed
- `fastcat` reverts to using a space separator (introduced in v0.16.0) between the Read ID and comment when outputting FASTQ comments that are not SAM tags

## [v0.16.7]
### Fixed
- Modification of BAM record with strtok when inferring Run ID from RG aux tag causing missing NM tag

## [v0.16.6]
### Fixed
- Additional spurious "contains non-integer 'NM' tag type" errors by checking EINVAL only when NM appears to be zero, and clearing errno first

## [v0.16.5]
### Fixed
- Spurious "contains non-integer 'NM' tag type" errors by checking EINVAL only when NM appears to be zero

## [v0.16.4]
### Added
- `bamstats` now saves histograms for unmapped reads when `--unmapped` is provided.

## [v0.16.3]
### Fixed
- Incorrect sanity check of NM.

## [v0.16.2]
### Fixed
- Prevent reads with implausible NM tag leading to illegal memory access in add_qual_count

## [v0.16.1]
### Changed
- Extended FASTQ SAM tag parsing to comment lines that include the RD tag (as well as RG).

## [v0.16.0]
### Added
- Support for reading SAM tags from FASTQ headers.
### Changed
- `fastcat` will output a tab between the Read ID and the SAM tags rather than a space to match samtools convention.
- `bamstats` uses `bam_get_tag_caseinsensitive` wrapper to get SAM tags with case insensitivity.
- `fastcat` and `bamstats` will infer a Run ID from the `RG` tag if `RD` is not available.
- Bumped version of htslib used to 1.19.
### Fixed
- Incorrectly capitalised ONT SAM tags are now output in lowercase by fastcat: `ch`, `rn`, `st`.

## [v0.15.2]
### Fixed
- Duplicated recipe name in Makefile.

### Added
- Section explaining `bamstats` output columns to README.

## [v0.15.1]
### Fixed
- Decimal precision of hisotgram outputs.

## [v0.15.0]
### Added
- Calculation of read length and quality  histograms to `fastcat` and `bamstats`.
- Calculation of alignment accuracy and alignment read coverage to `bamstats`.

## [v0.14.1]
### Fixed
- Missing compilation of conda aarch64 package

## [v0.14.0]
### Added
- `bamstats --duplex` option allows to count the number of duplex reads and 
  duplex-forming reads.

## [v0.13.2]
### Fixed
- Bug writing long reads to demultiplexed gzipped outputs.

## [v0.13.1]
### Fixed
- Bug writing `UINTMAX_MAX` for `min_length` and `nan` for `mean_quality` of a
  file in fastcat per-file stats if there were no reads in that file.

## [v0.13.0]
### Added
- Column with start time from MinKNOW header to `bamstats` output.

### Changed
- `bamstats` now prints `mean_quality`, `iden`, and `acc` values with 2 decimal
  places instead of 3 (the reason being that `fastcat` already uses 2 decimal
  places for `mean_quality` and more precision is unnecessary).

## [v0.12.0]
### Added
- Column with run ID from MinKNOW header to `fastcat` per-read stats and
  `bamstats` output.

## [v0.11.2]
### Changed
- Reverted the change of the default value of the `start_time` field to an empty
  string (it had been set to `"2000-01-01T00:00:00Z"` in v0.11.1).

## [v0.11.1]
### Fixed
- Bug in `fastcat` per-read summary stats.

## [v0.11.0]
### Changed
- Bamstats can now be run without a BAM index.
- `fastcat -H` now wraps all known header fields into SAM tags regardless of
  whether the header was "valid" (i.e. all expected fields were present) or not.

## [v0.10.2]
### Added
- Linux and macOS ARM conda packages.

## [v0.10.1]
### Fixed
- bamindex program missing from conda package.

## [v0.10.0]
### Added
- Create bamindex program to index unaligned BAMs for horizontal-parallel processing.

## [v0.9.0]
### Fixed
- Ensure reheadered fastq is indeed formatted as a valid SAM tag(s).

## [v0.8.0]
### Added
- Option to bamstats to add 'sample_name' column equivalent to fastcat.

## [v0.7.0]
### Added
- Option to report unmapped alignments in per read and summary files.

## [v0.6.1]
### Fixed
- Min read length in per-file statistics.

## [v0.6.0]
### Added
- `mean_quality` column to bamstats output, equivalent to that from fastcat.
- optional per-reference summary file for bamstats similar to samtools flagstats.

## [v0.5.0]
### Changed
- Behaviour of `-x/--recurse`. Top-level directory input will always be searched for
  data. Turning on recursion now exclusively refers to descending into child (and
  subsequent) directories.

## [v0.4.12]
### Fixed
- Updated kseq.h to allow exit on broken fastq/a stream.

## [v0.4.11]
### Changed
- `fastcat` will exit non-zero if an input file (named or recursed) cannot be opened

## [v0.4.10]
### Fixed
- Use of uninitialized memory in thread pool init, leading to memory leak.

## [v0.4.9]
### Fixed
- Handle BAM_CEQUAL and BAM_CDIFF that some aligners like to use.

## [v0.4.8]
### Fixed
- Doubled tab in output header.

## [v0.4.7]
### Changed
- Build conda package using bioconda's htslib.
### Fixed
- Occasional hanging on exit.

## [v0.4.6]
### Fixed
- Missing tab character in output header.

## [v0.4.5]
### Changed
- Pin openssl version in conda build to that which work across Python versions.

## [v0.4.4]
### Fixed
- Removed libdeflate from conda build which caused issues with threading.

## [v0.4.3]
### Changed
- Only multithread BAM decompression.

## [v0.4.2]
### Added
- Multithreading to `bamstats` for improved throughput.

## [v0.4.1]
### Changed
- Improved performance of `bamstats` for many-target bams.

## [v0.4.0]
### Added
- `bamstats` program for summarising (primary) alignment information.

## [v0.3.8]
### Fixed
- Refomatted header tags were space separated, fixed to tab separated.

## [v0.3.7]
### Added
- Option to reformat fastq headers as SAM-style tags for minimap2 passthrough.

## [v0.3.6]
### Fixed
- Per-file summary file created with broken header.

## [v0.3.5]
### Fixed
- Per-read summary file created incorrectly when `-s` option provided.

## [v0.3.4]
### Fixed
- Program hang when directory input given without trailing `/`.

## [v0.3.3]
### Added
- Transpose read number, channel, and start time from fastq headers to summary.
### Changed
- Additional columns in per-read summary file as above. These will be present,
  regardless of whether header information is present or not.

## [v0.3.2]
### Fixed
- Changed erroneously small MAX_BARCODE define; added runtime check to avoid
  invalid memory access.

## [v0.3.1]
### Fixed
- Updated CI release scripts.

## [v0.3.0]
### Added
- Parsing Guppy/MinKNOW fastq key=value header comments.
- Ability to demultiplex inputs based on "barcode" key in headers.
### Changed
- Per-read and per-file summary files now optional.

## [v0.2.1]
### Added
- Read length and read quality output filtering.
### Changed
- Average qualities computed with Kahan summation.

## [v0.2.0]
### Fixed
- Program hang when input file was non-existent or a directory.
### Added
- Ability to traverse a directory input.

## [v0.1.0]
### Added
- Ability to read input files from stdin.

## [v0.0.3]
### Changed
- Moved output files to optional arguments.
### Added
- `-s` option to add in a `sample_name` column to outputs.

## [v0.0.2]
### Changed
- No end-user changes.

## [v0.0.1]
### Added
- Per-read and per-file summarising of fastq data.

