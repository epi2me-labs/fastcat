# Contributing to fastcat2

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

- Use the `xalloc()`, `xrealloc()`, and `xrecalloc()` wrappers for allocations
- These wrappers log errors on allocation failure and exit immediately
- Free all allocated memory in cleanup/destroy functions
- For frequently-called hot-path code, consider reusable buffers to avoid per-call allocations

### Safety Checks

- Validate path lengths against `PATH_MAX`
- Use strict parsing functions (`fc_parse_*_strict()`) instead of `atoi()`/`atof()` for user input
- Clamp values to valid ranges when accessing lookup tables

### Function Scope

- Mark file-local functions as `static`
- Only expose functions in headers that are part of the public API
