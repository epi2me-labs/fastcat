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

- `output`: writes FASTQ or BAM and supports demultiplexed output.
- `perread`: writes per-read TSV; can include alignment columns when enabled.
- `perfile`: writes per-file TSV summary of read stats and failure counts.
- `counter`: writes run-id/basecaller counts per input file.
- `hist`: writes histogram tables from `fc_stats_engine` instances.
- `length_qscore_2d`: writes length/qscore heatmap tables from the writer model.
- `time_bin`: writes time-binned summaries from the writer model.
- `coverage`: forwards reads to the coverage writer.

## Example wiring

```c
fc_stats_engine_config stats_cfg = {
    .track_alignment = false,
    .track_polya = false,
    .polya_cover = 0.0f,
    .polya_qual = 0.0f,
    .polya_rev = false
};

fc_length_qscore_2d_config length_qscore_cfg = {
    .kll_k = 400,
    .qgrid_points = 257,
    .smoothing_points = 5
};

fc_time_bin_config time_bin_cfg = {
    .bin_size_s = 3600.0,
    .kll_k = 200
};

fc_writer_config writer_cfg = {
    .use_barcodes = false,
    .stats = &stats_cfg,
    .length_qscore = &length_qscore_cfg,
    .time_bin = &time_bin_cfg
};

fc_writer *writer = fc_writer_create(&writer_cfg);

fc_output_sink_opts out_opts = {
    .output_root = args.output_dir,   // NULL to stream FASTQ/BAM to stdout
    .reads_per_file = args.reads_per_file,
    .write_bam = args.write_bam,
    .threads = args.threads,
    .barcode_count = writer->barcode_count
};
fc_writer_add_sink(writer, fc_output_sink_create(&out_opts));

fc_perread_sink_opts perread_opts = {
    .output_root = args.output_dir,
    .sample = args.sample,
    .include_alignment = false,
    .use_barcodes = false,
    .barcode_count = 1
};
fc_writer_add_sink(writer, fc_perread_sink_create(&perread_opts));
```
