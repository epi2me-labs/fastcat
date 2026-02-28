# Contributing to fastcat

Thank you for your interest in contributing to fastcat!
This document outlines the development philosophy, coding standards, and contribution guidelines for this project.

## Development Philosophy

fastcat follows the **simple and transparent** design philosophy:
- Code should be straightforward and easy to understand
- Avoid unnecessary abstractions and over-engineering
- Prefer explicit code over "clever" solutions
- Only add complexity when directly required by the task at hand

## Error Handling Strategy

### Fail-Fast Approach

fastcat uses a **fail-fast** error handling strategy for convenience and clarity (and laziness):

- **Fatal errors**: When encountering unrecoverable errors or invalid input, programs should exit immediately with a clear error message
- Use `FC_LOG_FATAL(...)` for fatal errors that should terminate the program
- Use `exit(EXIT_FAILURE)` after logging critical errors
- This approach prevents silent failures and makes debugging easier in production workflows

### When to Fail Hard

Exit immediately for:
- Invalid or corrupted input data (e.g., missing required tags like NM for aligned reads)
- File I/O errors that prevent core functionality
- Memory allocation failures
- Configuration errors or invalid command-line arguments
- Violations of required preconditions (e.g., output directory already exists when fresh directory is required)

### Error Logging

Use the structured logging macros from `fc_log.h`:
- `FC_LOG_FATAL(...)` - Log fatal error and exit immediately
- `FC_LOG_ERROR(...)` - Log error message (does not exit)
- `FC_LOG_WARN(...)` - Log warning for recoverable issues
- `FC_LOG_INFO(...)` - Log informational messages
- `FC_LOG_DEBUG(...)` - Log debug information

Always provide actionable error messages that help users understand:
1. What went wrong
2. Why it's a problem
3. How to fix it (when applicable)

## Threading Model

fastcat is designed for **single-threaded record processing** (because the compute bottleneck is fastq/bam decompression):

- The main record processing pipeline (iterators, readers, writers, sinks) is NOT thread-safe
- Only the underlying BAM I/O layer (via htslib thread pool) may use multiple threads
- Do not add synchronization primitives (mutexes, atomics) for "future thread-safety"
- If adding code with lazy initialization or shared state, document that it requires single-threaded use

## Code Style

### General Guidelines

- Follow existing code style and patterns in the codebase
- Use meaningful variable and function names
- Keep functions focused and reasonably sized
- Add comments for non-obvious logic or complex algorithms

### Memory Management

- Use the `fc_alloc()`, `fc_realloc()`, and `fc_recalloc()` wrappers for allocations
- These wrappers log errors on allocation failure and exit immediately
- Free all allocated memory in cleanup/destroy functions
- For frequently-called hot-path code, consider reusable buffers to avoid per-call allocations

### String Formatting

Use project formatting helpers consistently. This keeps error handling uniform and avoids repeated manual buffer math.

- Do not call `snprintf()`, `asprintf()`, or `vasprintf()` directly in normal feature code
- For formatting into an existing/fixed buffer, use `FC_SNPRINTF_CHECK(...)`
- For dynamically allocated formatted strings, use `fc_xasprintf(...)` and free the result in normal cleanup paths
- Do not use `snprintf(NULL, 0, ...)` probe patterns at call sites; use `fc_xasprintf(...)` instead
- `vsnprintf(...)` is only acceptable in low-level varargs bridge code where a `va_list` is already required and performance rationale is explicit

Examples:

```c
char tmp[32];
FC_SNPRINTF_CHECK(tmp, sizeof(tmp), "barcode%04zu", barcode);

char *path = fc_xasprintf("output path", "%s/%s.%s", dir, stem, ext);
/* ... use path ... */
free(path);
```

Rationale:

- Centralized failure behavior (fail-fast with clear logging)
- Fewer manual length/overflow mistakes
- More readable call sites

### Safety Checks

- Validate path lengths against `PATH_MAX`
- Use strict parsing functions (`fc_parse_*_strict()`) instead of `atoi()`/`atof()` for user input
- Clamp values to valid ranges when accessing lookup tables

### Function Scope

- Mark file-local functions as `static`
- Only expose functions in headers that are part of the public API


## Code Audit

Instructions for automated code auditting tools

* Goal / window: Prepare the repo for an LTS-style release (~24 months) with minimal churn later; make necessary fixes/changes now to avoid risky/breaking work in the future.
* Scope / baseline: Review first-party, tracked files only on the dev branch; ignore untracked notes/TODOs; no open issues to incorporate.
* Review mode: Do a fresh, independent review each time (don’t anchor on prior findings).
* Focus area: Give extra scrutiny to FASTQ and BAM stats-generation subcommands.
* Priority order: Correctness → simplicity/readability → backward+forward compatibility → performance.
* Data model / security posture: Optimize for correctness on normal bioinformatics data (not primarily hostile/corrupt-file hardening).
* Input handling: Fail-fast on malformed input (fatal; no warn-and-continue).
* Compatibility: Treat CLI and outputs as part of the contract; aim for backward + forward compatibility over the LTS window. Expect no major FASTQ/BAM format changes, but allow for new meta tags that may later feed into summary stats.
* Platforms: Target Linux x86 via conda (per GitLab CI) for production; macOS is a primary dev environment.
* Readability vs performance: Prefer readability for non-C maintainers; do not undo intentional performance optimizations—document them if they’re non-obvious.
* Testing: Prefer unit tests; Makefile integration tests are acceptable but aim to reduce/curtail them where practical.
* Reproducibility: Tighten/pin dependencies/tooling to support reproducible builds.
* Change discipline: Prefer minimal, local, defensible changes; avoid speculative overengineering; make proactive fixes where justified.
* Numeric robustness: Handle “janky” values safely at the point of use (e.g., histogram binning, accumulation) without broad unlikely-case machinery unless justified.

Deliverables
* Output a comprehensive but ordered report as an issue list first, ordered by severity:
	1.	Critical (security/crashes/data corruption)
	2.	High (correctness bugs/edge cases)
	3.	Medium (maintainability/simplification)
	4.	Low (style/docs/optimization)
* Each issue must include: file/line, impact, and a concrete fix.
* Put questions/assumptions second and summary last.
* Keep commits small, test-backed, and independently reviewable.
* Run the repo’s static analysis and memory checks as part of the audit gate.
