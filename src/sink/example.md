# Sinks overview

Sinks are small writer backends that receive fully parsed reads and emit output.
`fc_writer` owns sink instances and calls these hooks:

- `begin_file(ctx, fname)`
- `add_read(ctx, read, read_ctx, fname, hdr)`
- `end_file(ctx, fname, stats)`
- `finish(ctx)`
- `destroy(ctx)`

Sinks should not implement demultiplexing or alignment policy. Those decisions are
made in `fc_writer` and encoded in `fc_read_ctx` for each read.

## Included sinks

- `output`: writes FASTQ or BAM, supports demultiplexed output and reheadering.
- `perread`: writes per-read TSV; can include alignment columns when enabled.
- `perfile`: writes per-file TSV summary of read stats and failure counts.
- `counter`: writes run-id/basecaller counts per input file.
- `hist`: writes histogram tables from `fc_stats_engine` instances.
- `coverage`: forwards reads to the coverage writer.

## Example wiring

```c
fc_writer *writer = fc_writer_create();

fc_output_sink_opts out_opts = {
    .output_dir = args.demultiplex_dir,
    .reads_per_file = args.reads_per_file,
    .reheader = args.reheader,
    .write_bam = args.write_bam,
    .threads = args.threads
};
fc_writer_add_sink(writer, fc_output_sink_create(&out_opts));

fc_perread_sink_opts perread_opts = {
    .path = args.perread,
    .sample = args.sample,
    .include_alignment = false,
    .output_unmapped = true
};
fc_writer_add_sink(writer, fc_perread_sink_create(&perread_opts));
```
