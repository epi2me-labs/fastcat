# fastcat

A simple utility to concatenate .fastq(.gz) files whilst creating a summary
of the sequences:

```
./fastcat output.txt reads1.fastq(.gz) reads2.fastq(.gz) ... | gzip > all_reads.fastq.gz
```

The `output.txt` is a tab-seaparate file with columns:

```
read_id        filename                read_length  mean_quality
SRR12447496.1  SRR12447496_1.fastq.gz  531          14.03
SRR12447496.2  SRR12447496_1.fastq.gz  513          13.91
SRR12447496.3  SRR12447496_1.fastq.gz  473          14.70
...
```

The mean quality is defined as:
```
-10 * log10(mean(10^(Q/-10)))
```

where `Q` are the set of all per-base quality scores for the read.
