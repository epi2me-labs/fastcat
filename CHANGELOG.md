# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


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

