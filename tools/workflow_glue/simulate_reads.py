"""Stream synthetic FASTQ or aligned-SAM records to stdout.

Designed for fastcat testing:
- timestamps sampled from a throughput curve over time
- read classes sampled from a small latent mixture
- read length, qscore, and translocation speed sampled jointly by class
- explicit SAM-style `qs:f:` tag in FASTQ comment so fastcat can use it
  directly without recomputing quality
- optional aligned-SAM output with varied primary, secondary, supplementary,
  duplicate, QC-fail, and unmapped records for report demonstrations

Example:
  python3 tools/dev-tools simulate_reads \
    --n 200000 --shape 5.0 --scale 6000 \
    --q1 11 --q2 28 \
    --barcodes barcode0001 barcode0002 unclassified \
    --barcode-weights 0.7 0.2 0.1 \
    > /tmp/synth.fastq
"""

import argparse
import bisect
import csv
import datetime as dt
import math
from pathlib import Path
import random
import sys

import pysam


DEFAULT_RUNID = "12345678-1234-1234-1234-123456789abc"
DEFAULT_MODEL = "dna_r10.4.1_e8.2_400bps_hac@v5.0.0"
DEFAULT_FLOW_CELL_ID = "XXX12345"
DEFAULT_SPEED_MEAN = 400.0
DEFAULT_SPEED_IQR = 40.0
DEFAULT_SPEED_FINAL_MEAN = 340.0
DEFAULT_THROUGHPUT_BIN_S = 1800.0
SPEED_QSCORE_BAND_MIN = 390.0
SPEED_QSCORE_BAND_MAX = 420.0
SPEED_QSCORE_MAX_PENALTY = 5.0
SPEED_QSCORE_PENALTY_SPAN = 60.0
CHANNEL_MAP_ROOT = Path(__file__).resolve().parent / "data"
CHANNEL_MAP_PATHS = {
    "flongle": CHANNEL_MAP_ROOT / "asicb512fx_flongle_channel_map.csv",
    "minion": CHANNEL_MAP_ROOT / "asicb512fx_channel_map.csv",
    "promethion": CHANNEL_MAP_ROOT / "3000_opx_channel_map.csv",
}
DEVICE_CHANNELS = {
    "flongle": 128,
    "minion": 512,
    "promethion": 3000,
}
READ_CLASSES = ("good", "mediocre", "junk_short")


class RunningStats:
    def __init__(self):
        self.count = 0
        self.sum_len = 0.0
        self.sum_q = 0.0
        self.min_len = 1 << 62
        self.max_len = 0

    def add(self, length, q):
        self.count += 1
        self.sum_len += float(length)
        self.sum_q += float(q)
        if length < self.min_len:
            self.min_len = length
        if length > self.max_len:
            self.max_len = length


def draw_from_mix(rng, mix):
    r = rng.random()
    c = 0.0
    for name, p in mix:
        c += p
        if r <= c:
            return name
    return mix[-1][0]


def clamp(x, lo, hi):
    if x < lo:
        return lo
    if x > hi:
        return hi
    return x


def smoothstep(x):
    x = clamp(x, 0.0, 1.0)
    return x * x * (3.0 - 2.0 * x)


def gaussian_bump(x, centre, width):
    if width <= 0.0:
        return 0.0
    z = (x - centre) / width
    return math.exp(-0.5 * z * z)


def throughput_profile(frac):
    # quick startup, long decay, two outages, slight recovery late
    if frac < 0.01:
        rate = 0.18 + 75.0 * frac
    elif frac < 0.05:
        rate = 0.93 + 0.04 * math.sin(18.0 * math.pi * frac)
    elif frac < 0.22:
        rate = 0.95 - 0.25 * ((frac - 0.05) / 0.17)
    elif frac < 0.55:
        rate = 0.77 - 0.18 * ((frac - 0.22) / 0.33)
    elif frac < 0.75:
        rate = 0.72 - 0.12 * ((frac - 0.55) / 0.20)
    else:
        tail = (frac - 0.75) / 0.25
        rate = 0.60 * math.exp(-3.6 * tail) + 0.012

    rate *= 1.0 - 0.98 * gaussian_bump(frac, 0.245, 0.009)
    rate *= 1.0 - 0.995 * gaussian_bump(frac, 0.455, 0.008)
    rate *= 1.0 + 0.18 * gaussian_bump(frac, 0.83, 0.035)
    return max(rate, 1e-6)


def build_time_sampler(args):
    duration = timestamp_duration_seconds(args)
    if args.time_profile == "linear" or duration <= 0.0:
        return None

    bins = max(1, int(math.ceil(duration / DEFAULT_THROUGHPUT_BIN_S)))
    cumulative = []
    total = 0.0
    for idx in range(bins):
        mid_s = min(duration, (idx + 0.5) * DEFAULT_THROUGHPUT_BIN_S)
        frac = mid_s / duration if duration > 0.0 else 0.0
        total += throughput_profile(frac)
        cumulative.append(total)
    return cumulative, total, duration


def sample_timestamp_from_profile(rng, sampler):
    cumulative, total, duration = sampler
    if total <= 0.0:
        return 0.0
    target = rng.random() * total
    idx = bisect.bisect_left(cumulative, target)
    if idx >= len(cumulative):
        idx = len(cumulative) - 1
    start = idx * DEFAULT_THROUGHPUT_BIN_S
    end = min(duration, start + DEFAULT_THROUGHPUT_BIN_S)
    if end <= start:
        return start
    return start + rng.random() * (end - start)


def class_mix_for_time(frac):
    late = smoothstep((frac - 0.68) / 0.32)
    outage = max(
        gaussian_bump(frac, 0.245, 0.018),
        gaussian_bump(frac, 0.455, 0.016),
    )
    recovery = gaussian_bump(frac, 0.83, 0.06)

    good = 0.84 - 0.14 * late - 0.06 * outage + 0.03 * recovery
    mediocre = 0.12 + 0.08 * late + 0.03 * outage
    junk = 1.0 - good - mediocre
    junk = clamp(junk, 0.008, 0.06)
    total = good + mediocre + junk
    return [
        ("good", good / total),
        ("mediocre", mediocre / total),
        ("junk_short", junk / total),
    ]


def good_speed_mean(start_mean, final_mean, frac):
    if frac < 0.70:
        base = start_mean + 8.0 * smoothstep(frac / 0.70)
    else:
        late = (frac - 0.70) / 0.30
        top = start_mean + 8.0
        base = top + (final_mean - top) * smoothstep(late)
    base -= 10.0 * gaussian_bump(frac, 0.245, 0.022)
    base -= 8.0 * gaussian_bump(frac, 0.455, 0.020)
    base += 8.0 * gaussian_bump(frac, 0.83, 0.05)
    return base


def time_quality_offset(frac):
    late = smoothstep((frac - 0.72) / 0.28)
    return -1.4 * late - 1.1 * gaussian_bump(frac, 0.245, 0.020) - 0.9 * gaussian_bump(frac, 0.455, 0.018)


def time_noise_scale(frac):
    return 1.0 + 0.85 * smoothstep((frac - 0.72) / 0.28)


def class_length(rng, read_class, args):
    if read_class == "good":
        if rng.random() < 0.84:
            length = rng.gammavariate(args.shape, args.scale)
        else:
            length = rng.gammavariate(2.0, 1600.0)
    elif read_class == "mediocre":
        if rng.random() < 0.65:
            length = rng.gammavariate(3.4, 2200.0)
        else:
            length = rng.gammavariate(1.8, 5200.0)
    else:
        if rng.random() < 0.82:
            length = rng.gammavariate(1.25, 260.0)
        else:
            length = rng.gammavariate(1.3, 1400.0)
    return int(round(clamp(length, args.min_len, args.max_len)))


def class_speed(rng, read_class, args, frac):
    good_mean = good_speed_mean(args.speed_mean, args.speed_final_mean, frac)
    sd_scale = time_noise_scale(frac)
    if read_class == "good":
        mean = good_mean
        sd = args.speed_sd * (0.55 + 0.25 * sd_scale)
    elif read_class == "mediocre":
        mean = good_mean - (22.0 + 10.0 * smoothstep((frac - 0.45) / 0.55))
        sd = args.speed_sd * (0.90 + 0.25 * sd_scale)
    else:
        mean = 275.0 - 35.0 * smoothstep((frac - 0.60) / 0.40)
        sd = 12.0 + 8.0 * sd_scale
    return sample_speed(rng, mean, sd)


def class_qscore(rng, read_class, length, speed_nt_s, args, frac):
    length_kb = float(length) / 1000.0
    speed_penalty = qscore_speed_penalty(speed_nt_s)
    time_shift = time_quality_offset(frac)
    noise_scale = time_noise_scale(frac)

    if read_class == "good":
        base = args.q2 + 0.3 + time_shift
        mixture = rng.random()
        if mixture < 0.44:
            base -= rng.uniform(1.0, 8.0)
        elif mixture > 0.94:
            base += rng.uniform(6.5, 9.0)
        sigma = args.q2_sd * (0.72 + 0.18 * noise_scale)
        q = rng.gauss(base - 0.45 * speed_penalty, sigma)
    elif read_class == "mediocre":
        base = args.q1 + 0.4 + 0.30 * time_shift
        base += 0.10 * min(length_kb, 12.0)
        if rng.random() < 0.25:
            base += rng.uniform(1.0, 3.0)
        sigma = args.q1_sd * (0.85 + 0.20 * noise_scale)
        q = rng.gauss(base - 0.60 * speed_penalty, sigma)
    else:
        base = 7.2 - 0.4 * smoothstep((frac - 0.65) / 0.35)
        sigma = 0.75 + 0.15 * noise_scale
        q = rng.gauss(base - 0.18 * speed_penalty, sigma)

    return clamp(q, args.q_min, args.q_max)


def load_channel_centres(device, channels):
    centres = []
    path = CHANNEL_MAP_PATHS[device]
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, skipinitialspace=True)
        for row in reader:
            channel = int(row["name"])
            if channel > channels:
                continue
            xs = []
            ys = []
            for mux in range(1, 5):
                x_key = f"mux_{mux}_phys_x"
                y_key = f"mux_{mux}_phys_y"
                if x_key not in row or y_key not in row:
                    x_key = "mux_1_phys_x"
                    y_key = "mux_1_phys_y"
                xs.append(float(row[x_key]))
                ys.append(float(row[y_key]))
            x = sum(xs) / len(xs)
            y = sum(ys) / len(ys)
            if device == "minion":
                x, y = y, x
            centres.append((channel, x, y))
    if not centres:
        centres.append((1, 0.0, 0.0))
    return centres


def build_channel_sampler(device, channels):
    centres = load_channel_centres(device, channels)
    xs = [x for _, x, _ in centres]
    ys = [y for _, _, y in centres]
    x_mid = 0.5 * (min(xs) + max(xs))
    y_mid = 0.5 * (min(ys) + max(ys))
    x_radius = max(1e-6, 0.46 * (max(xs) - min(xs)))
    y_radius = max(1e-6, 0.46 * (max(ys) - min(ys)))

    cumulative = []
    channel_ids = []
    total = 0.0
    for channel, x, y in centres:
        dx = (x - x_mid) / x_radius
        dy = (y - y_mid) / y_radius
        r2 = dx * dx + dy * dy
        weight = 1.0 + 0.35 * math.exp(-0.75 * r2)
        total += weight
        cumulative.append(total)
        channel_ids.append(channel)
    return channel_ids, cumulative, total


def sample_channel(rng, channel_sampler):
    channel_ids, cumulative, total = channel_sampler
    target = rng.random() * total
    idx = bisect.bisect_left(cumulative, target)
    if idx >= len(channel_ids):
        idx = len(channel_ids) - 1
    return channel_ids[idx]


def phred_char(q):
    # FASTQ Sanger: ASCII 33..126 (practically we clamp phred 0..60)
    qi = int(round(clamp(q, 0.0, 60.0)))
    return chr(33 + qi)


def make_seq(rng, length, random_seq):
    if not random_seq:
        return "A" * length
    bases = "ACGT"
    return "".join(bases[rng.getrandbits(2)] for _ in range(length))


def sample_speed(rng, mean, sd):
    return max(1e-3, rng.gauss(mean, sd))


def qscore_speed_penalty(speed_nt_s):
    if SPEED_QSCORE_BAND_MIN <= speed_nt_s <= SPEED_QSCORE_BAND_MAX:
        return 0.0
    if speed_nt_s < SPEED_QSCORE_BAND_MIN:
        distance = SPEED_QSCORE_BAND_MIN - speed_nt_s
    else:
        distance = speed_nt_s - SPEED_QSCORE_BAND_MAX
    penalty = SPEED_QSCORE_MAX_PENALTY * (distance / SPEED_QSCORE_PENALTY_SPAN) ** 2
    return min(SPEED_QSCORE_MAX_PENALTY, penalty)


def iso8601_with_seconds(base, seconds):
    return (base + dt.timedelta(seconds=seconds)).strftime("%Y-%m-%dT%H:%M:%SZ")


def timestamp_duration_seconds(args):
    if args.run_duration_s is not None:
        return args.run_duration_s
    if args.n <= 1:
        return 0.0
    return float(args.n - 1) * float(args.seconds_step)


def sample_timestamp_offset_seconds(rng, args, read_index):
    if args.time_profile == "linear":
        if read_index <= 1:
            return 0.0
        return float(read_index - 1) * float(args.seconds_step)

    duration = timestamp_duration_seconds(args)
    if duration <= 0.0:
        return 0.0

    u = rng.random()
    tau = duration * args.time_tau_fraction
    cutoff = 1.0 - math.exp(-duration / tau)
    return -tau * math.log1p(-u * cutoff)


def alignment_error_count(rng, aligned_length, mean_q, args):
    """Sample an NM count from the simulated read quality."""
    concordance_qscore = max(0.0, mean_q)
    error_probability = 10.0 ** (-concordance_qscore / 10.0)
    expected = aligned_length * error_probability
    standard_deviation = math.sqrt(
        aligned_length * error_probability * (1.0 - error_probability))
    return max(0, min(
        aligned_length - 1,
        int(round(rng.gauss(expected, standard_deviation)))))


def primary_alignment_cigar(rng, length, mean_q, args):
    """Create a valid clipped CIGAR and quality-calibrated NM tag."""
    if args.aligned_fraction_sd == 0.0:
        aligned_fraction = args.aligned_fraction_mean
    else:
        variance = args.aligned_fraction_sd ** 2
        concentration = (
            args.aligned_fraction_mean
            * (1.0 - args.aligned_fraction_mean)
            / variance
            - 1.0)
        aligned_fraction = rng.betavariate(
            args.aligned_fraction_mean * concentration,
            (1.0 - args.aligned_fraction_mean) * concentration)
    aligned_fraction = max(args.aligned_fraction_min, aligned_fraction)
    clipped = int(round(length * (1.0 - aligned_fraction)))
    left_clip = rng.randint(0, clipped)
    right_clip = clipped - left_clip
    query_aligned = length - clipped

    errors = alignment_error_count(rng, query_aligned, mean_q, args)
    max_insertion = min(
        errors // 4,
        int(query_aligned * args.max_indel_fraction),
        int(query_aligned * args.insertion_rate_max))
    insertion = rng.randint(0, max_insertion)
    match = query_aligned - insertion
    max_deletion = min(
        (errors - insertion) // 4,
        int(match * args.max_indel_fraction),
        int(match * args.deletion_rate_max))
    deletion = rng.randint(0, max_deletion)
    split = rng.randint(max(1, match // 4), max(1, 3 * match // 4))
    parts = []
    if left_clip:
        parts.append(f"{left_clip}S")
    parts.append(f"{split}M")
    if insertion:
        parts.append(f"{insertion}I")
    parts.append(f"{match - split}M")
    if deletion:
        parts.append(f"{deletion}D")
    if right_clip:
        parts.append(f"{right_clip}S")

    return "".join(parts), errors, match + deletion


def sam_header(args):
    """Build the compact two-reference header used for simulated alignments."""
    read_groups = []
    for index, barcode in enumerate(args.barcodes):
        group = {
            "ID": f"simulated_{index}",
            "DS": f"runid={args.runid} basecall_model={args.model}",
            "LB": args.sample_id or "simulated-data",
            "PL": "ONT",
            "PU": args.flow_cell_id,
        }
        if barcode != "unclassified":
            group["SM"] = barcode
        read_groups.append(group)
    return pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate" if args.sam_coordinate_sort else "unsorted"},
        "SQ": [
            {"SN": "chr1", "LN": args.sam_reference_length},
            {"SN": "chr2", "LN": args.sam_reference_length},
        ],
        "RG": read_groups,
    })


def read_group_id(args, barcode):
    return f"simulated_{args.barcodes.index(barcode)}"


def fastq_read_group_id(args, barcode):
    suffix = "" if barcode == "unclassified" else f"_{barcode}"
    return f"{args.runid}_{args.model}{suffix}"


def sam_record(header, args, name, flag, reference, position, cigar, sequence, quality, tags):
    """Build one aligned or unaligned record through pysam."""
    record = pysam.AlignedSegment(header)
    record.query_name = name
    record.flag = flag
    record.query_sequence = sequence
    record.query_qualities = pysam.qualitystring_to_array(quality)
    record.set_tags(tags)
    if reference == "*":
        record.reference_id = -1
        record.reference_start = -1
        record.mapping_quality = 0
        record.cigarstring = None
    else:
        record.reference_id = header.get_tid(reference)
        record.reference_start = position - 1
        record.mapping_quality = args.sam_mapping_quality
        record.cigarstring = cigar
    return record


def emit_sam_record(writer, records, record):
    """Write one SAM record now, or retain it for coordinate sorting."""
    if records is None:
        writer.write(record)
        return
    records.append(record)


def write_fastq_record(out, record, args, barcode, omit_metadata):
    """Write an unaligned pysam record using ONT's FASTQ header convention."""
    if omit_metadata:
        out.write(f"@{record.query_name}\n")
    else:
        def fastq_tag(name, value, kind):
            # pysam uses the narrowest BAM integer representation (for example
            # C for small values); FASTQ comments use the portable SAM i type.
            if kind in "cCsSiI":
                kind = "i"
            return f"{name}:{kind}:{value}"

        record_tags = [
            fastq_tag(name, value, kind)
            for name, value, kind in record.get_tags(with_value_type=True)
            if name != "RG"
        ]
        header_tags = [
            f"RG:Z:{fastq_read_group_id(args, barcode)}",
            f"DT:Z:{args.start_time}",
            f"PU:Z:{args.flow_cell_id}",
        ]
        if args.sample_id:
            header_tags.append(f"LB:Z:{args.sample_id}")
        if barcode != "unclassified":
            header_tags.append(f"SM:Z:{barcode}")
        tags = "\t".join(header_tags + record_tags)
        out.write(f"@{record.query_name} {tags}\n")
    out.write(f"{record.query_sequence}\n+\n{record.qual}\n")


def argparser():
    p = argparse.ArgumentParser(
        description="Stream synthetic FASTQ with gamma lengths and bimodal mean-q; emits qs:f tag.",
        add_help=False,
    )

    run = p.add_argument_group("simulation")
    run.add_argument("--n", type=int, default=10000, help="Number of reads to emit (default: 10000)")
    run.add_argument("--seed", type=int, default=42, help="RNG seed (default: 42)")

    length = p.add_argument_group("read length distribution")
    length.add_argument("--shape", type=float, default=5.0, help="Main productive-read Gamma shape k (default: 5.0)")
    length.add_argument("--scale", type=float, default=6000.0, help="Main productive-read Gamma scale theta (default: 6000)")
    length.add_argument("--min-len", type=int, default=100, help="Minimum read length clamp (default: 100)")
    length.add_argument("--max-len", type=int, default=1_000_000, help="Maximum read length clamp (default: 1000000)")

    qscore = p.add_argument_group("qscore distribution")
    qscore.add_argument("--q1", type=float, default=12.5, help="Mediocre-read qscore mean (default: 12.5)")
    qscore.add_argument("--q2", type=float, default=20.5, help="Good-read qscore mean (default: 20.5)")
    qscore.add_argument("--q1-sd", type=float, default=2.1, help="Stddev around mediocre-read qscore mean (default: 2.1)")
    qscore.add_argument("--q2-sd", type=float, default=2.6, help="Stddev around good-read qscore mean (default: 2.6)")
    qscore.add_argument("--q1-weight", type=float, default=0.15, help="Legacy option; retained for compatibility (default: 0.15)")
    qscore.add_argument("--q-min", type=float, default=1.0, help="Minimum qscore clamp (default: 1)")
    qscore.add_argument("--q-max", type=float, default=40.0, help="Maximum qscore clamp (default: 40)")

    barcode = p.add_argument_group("barcode mixture")
    barcode.add_argument("--barcodes", nargs="+", default=["unclassified"], help="Barcode names to sample")
    barcode.add_argument("--barcode-weights", nargs="+", type=float, default=[1.0], help="Weights matching --barcodes")

    metadata = p.add_argument_group("read metadata")
    metadata.add_argument("--runid", default=DEFAULT_RUNID, help=f"Run ID for RG composition (default: {DEFAULT_RUNID})")
    metadata.add_argument("--model", default=DEFAULT_MODEL, help=f"Basecaller model in RG (default: {DEFAULT_MODEL})")
    metadata.add_argument(
        "--flow-cell-id", default=DEFAULT_FLOW_CELL_ID,
        help=f"Flow cell ID metadata to emit (default: {DEFAULT_FLOW_CELL_ID})")
    metadata.add_argument("--id-prefix", default="simread", help="Read ID prefix (default: simread)")
    metadata.add_argument("--sample-id", default="", help="Sample ID metadata to emit (default: unset)")
    metadata.add_argument(
        "--omit-metadata",
        action="store_true",
        help="Emit bare FASTQ headers without optional metadata tags")
    metadata.add_argument(
        "--device", choices=["flongle", "minion", "promethion"], default="minion",
        help="Flowcell layout to simulate (default: minion)")
    metadata.add_argument(
        "--channels", type=int, default=None,
        help="Channel range 1..N; defaults to the full range for --device")
    metadata.add_argument("--muxes", type=int, default=4, help="Mux range 1..N (default: 4)")
    metadata.add_argument("--start-time", default="2026-01-01T00:00:00Z", help="Base start_time tag (default: 2026-01-01T00:00:00Z)")
    metadata.add_argument(
        "--time-profile", choices=["saturating", "linear"], default="saturating",
        help="Timestamp spacing profile; 'saturating' now uses a realistic throughput curve (default: saturating)")
    metadata.add_argument(
        "--seconds-step", type=int, default=1,
        help="Linear timestamp step, and default run-duration basis for saturating profile (default: 1)")
    metadata.add_argument(
        "--run-duration-s", type=float, default=None,
        help="Final timestamp offset for saturating profile; defaults to (n - 1) * seconds-step")
    metadata.add_argument(
        "--time-tau-fraction", type=float, default=0.2,
        help="Truncated-exponential timestamp time constant as a fraction of run duration (default: 0.2)")
    metadata.add_argument("--duplex-rate", type=float, default=0.0, help="Fraction tagged as duplex via dx:i:1 (default: 0)")

    speed = p.add_argument_group("translocation speed")
    speed.add_argument(
        "--speed-mean", type=float, default=DEFAULT_SPEED_MEAN,
        help=f"Initial mean translocation speed in nt/s for du:f generation (default: {DEFAULT_SPEED_MEAN})")
    speed.add_argument(
        "--speed-final-mean", type=float, default=DEFAULT_SPEED_FINAL_MEAN,
        help=f"Final mean translocation speed in nt/s at the end of the run (default: {DEFAULT_SPEED_FINAL_MEAN})")
    speed.add_argument(
        "--speed-iqr", type=float, default=DEFAULT_SPEED_IQR,
        help=f"Interquartile spread of translocation speed in nt/s (default: {DEFAULT_SPEED_IQR})")
    speed.add_argument(
        "--speed-sd", type=float, default=None,
        help="Optional stddev of translocation speed in nt/s; overrides --speed-iqr if set")

    output = p.add_argument_group("output")
    output.add_argument(
        "--output-format", choices=["fastq", "sam"], default="fastq",
        help="Record format to write; use samtools to convert SAM to BAM (default: fastq)")
    output.add_argument("--random-seq", action="store_true", help="Emit random A/C/G/T instead of homopolymer A")
    output.add_argument("--progress-every", type=int, default=0, help="If >0, print progress to stderr every N reads")

    alignment = p.add_argument_group("SAM alignment simulation")
    alignment.add_argument(
        "--aligned-rate", type=float, default=0.99,
        help="Fraction of simulated reads with a primary alignment in SAM output (default: 0.99)")
    alignment.add_argument(
        "--secondary-rate", type=float, default=0.05,
        help="Fraction of primary alignments with an additional secondary record (default: 0.05)")
    alignment.add_argument(
        "--supplementary-rate", type=float, default=0.04,
        help="Fraction of primary alignments with an additional supplementary record (default: 0.04)")
    alignment.add_argument(
        "--reverse-rate", type=float, default=0.5,
        help="Fraction of primary alignments placed on the reverse strand (default: 0.5)")
    alignment.add_argument(
        "--qcfail-rate", type=float, default=0.025,
        help="Fraction of primary alignments flagged QC-fail (default: 0.025)")
    alignment.add_argument(
        "--duplicate-rate", type=float, default=0.04,
        help="Fraction of primary alignments flagged duplicate (default: 0.04)")
    alignment.add_argument(
        "--duplex-forming-rate", type=float, default=0.03,
        help="Fraction of non-duplex reads tagged duplex-forming (default: 0.03)")
    alignment.add_argument(
        "--second-reference-rate", type=float, default=1.0 / 3.0,
        help="Fraction of primary alignments placed on chr2 (default: 0.333)")
    alignment.add_argument(
        "--aligned-fraction-mean", type=float, default=0.956,
        help="Mean fraction of each mapped read represented by its alignment (default: 0.956)")
    alignment.add_argument(
        "--aligned-fraction-sd", type=float, default=0.06,
        help="Standard deviation of mapped-read alignment fraction (default: 0.06)")
    alignment.add_argument(
        "--aligned-fraction-min", type=float, default=0.55,
        help="Minimum mapped-read alignment fraction after clamping (default: 0.55)")
    alignment.add_argument(
        "--insertion-rate-max", type=float, default=0.025,
        help="Maximum simulated insertion fraction within a mapped read (default: 0.025)")
    alignment.add_argument(
        "--deletion-rate-max", type=float, default=0.02,
        help="Maximum simulated deletion fraction within a mapped read (default: 0.02)")
    alignment.add_argument(
        "--max-indel-fraction", type=float, default=1.0 / 12.0,
        help="Hard maximum fraction of aligned bases represented by one indel type (default: 0.0833)")
    alignment.add_argument(
        "--sam-mapping-quality", type=int, default=60,
        help="Mapping quality for simulated mapped SAM records (default: 60)")
    alignment.add_argument(
        "--sam-reference-length", type=int, default=5_000_000,
        help="Length of each synthetic SAM reference (default: 5000000)")
    alignment.add_argument(
        "--sam-coordinate-sort", action="store_true",
        help="Buffer SAM records and emit coordinate order; intended for small coverage demos")

    return p


def validate_args(p, args):
    if args.n < 0:
        p.error("--n must be >= 0")
    if args.shape <= 0.0 or args.scale <= 0.0:
        p.error("--shape and --scale must be > 0")
    if args.min_len <= 0 or args.max_len < args.min_len:
        p.error("--min-len must be > 0 and --max-len must be >= --min-len")
    if not (0.0 <= args.q1_weight <= 1.0):
        p.error("--q1-weight must be in [0,1]")
    if args.q1_sd <= 0.0 or args.q2_sd <= 0.0:
        p.error("--q1-sd and --q2-sd must be > 0")
    if args.q_min > args.q_max:
        p.error("--q-min must be <= --q-max")
    if args.channels is None:
        args.channels = DEVICE_CHANNELS[args.device]
    if args.channels <= 0:
        p.error("--channels must be > 0")
    if args.channels > DEVICE_CHANNELS[args.device]:
        p.error("--channels must not exceed the selected device capacity")
    if args.muxes <= 0:
        p.error("--muxes must be > 0")
    if args.seconds_step < 0:
        p.error("--seconds-step must be >= 0")
    if args.run_duration_s is not None and (
        not math.isfinite(args.run_duration_s) or args.run_duration_s < 0.0
    ):
        p.error("--run-duration-s must be finite and >= 0")
    if not math.isfinite(args.time_tau_fraction) or args.time_tau_fraction <= 0.0:
        p.error("--time-tau-fraction must be finite and > 0")
    if not (0.0 <= args.duplex_rate <= 1.0):
        p.error("--duplex-rate must be in [0,1]")
    if not math.isfinite(args.speed_mean) or args.speed_mean <= 0.0:
        p.error("--speed-mean must be finite and > 0")
    if not math.isfinite(args.speed_final_mean) or args.speed_final_mean <= 0.0:
        p.error("--speed-final-mean must be finite and > 0")
    if not math.isfinite(args.speed_iqr) or args.speed_iqr <= 0.0:
        p.error("--speed-iqr must be finite and > 0")
    if args.speed_sd is not None and (not math.isfinite(args.speed_sd) or args.speed_sd <= 0.0):
        p.error("--speed-sd must be finite and > 0")
    if args.progress_every < 0:
        p.error("--progress-every must be >= 0")
    for option in (
            "aligned_rate", "secondary_rate", "supplementary_rate", "reverse_rate",
            "qcfail_rate", "duplicate_rate", "duplex_forming_rate", "second_reference_rate",
            "aligned_fraction_mean", "aligned_fraction_min", "insertion_rate_max",
            "deletion_rate_max", "max_indel_fraction"):
        if not 0.0 <= getattr(args, option) <= 1.0:
            p.error(f"--{option.replace('_', '-')} must be in [0,1]")
    if args.aligned_fraction_sd < 0.0:
        p.error("--aligned-fraction-sd must be >= 0")
    if (
            args.aligned_fraction_sd > 0.0
            and args.aligned_fraction_sd ** 2
            >= args.aligned_fraction_mean * (1.0 - args.aligned_fraction_mean)):
        p.error("--aligned-fraction-sd is too large for --aligned-fraction-mean")
    if args.aligned_fraction_min > 1.0:
        p.error("--aligned-fraction-min must be <= 1")
    if not 0 <= args.sam_mapping_quality <= 255:
        p.error("--sam-mapping-quality must be in [0,255]")
    if (
            args.output_format == "sam"
            and args.sam_reference_length < math.ceil(
                args.max_len * (1.0 + args.deletion_rate_max))):
        p.error("--sam-reference-length must allow for --max-len plus simulated deletions")

    if len(args.barcodes) != len(args.barcode_weights):
        p.error("--barcodes and --barcode-weights must have the same number of entries")
    if any(not name.strip() for name in args.barcodes):
        p.error("--barcodes must not contain empty names")
    if any(not math.isfinite(w) or w < 0.0 for w in args.barcode_weights):
        p.error("--barcode-weights must be finite and >= 0")
    barcode_weight_total = sum(args.barcode_weights)
    if barcode_weight_total <= 0.0:
        p.error("--barcode-weights must sum to > 0")
    args.barcode_mix = [
        (name.strip(), weight / barcode_weight_total)
        for name, weight in zip(args.barcodes, args.barcode_weights)
    ]

    try:
        args.start_dt = dt.datetime.strptime(args.start_time, "%Y-%m-%dT%H:%M:%SZ")
    except ValueError:
        p.error("--start-time must be UTC ISO8601 like 2026-01-01T00:00:00Z")

    if args.speed_sd is None:
        args.speed_sd = args.speed_iqr / 1.3489795003921634


def main(args):
    validate_args(argparser(), args)
    rng = random.Random(args.seed)
    time_rng = random.Random(args.seed + 1_000_003)
    channel_rng = random.Random(args.seed + 2_000_003)
    channel_sampler = build_channel_sampler(args.device, args.channels)

    err = sys.stderr
    sam_header_value = sam_header(args)
    sam_writer = None
    sam_records = None
    if args.output_format == "sam":
        sam_writer = pysam.AlignmentFile("-", "w", header=sam_header_value)
        if args.sam_coordinate_sort:
            sam_records = []

    stats = RunningStats()
    class_counts = {name: 0 for name in READ_CLASSES}
    duplex_n = 0
    max_time_s = 0.0
    sum_speed = 0.0
    barcode_counts = {name: 0 for name, _ in args.barcode_mix}
    alignment_counts = {
        "primary": 0,
        "unmapped": 0,
        "secondary": 0,
        "supplementary": 0,
    }
    run_duration_s = timestamp_duration_seconds(args)
    time_sampler = build_time_sampler(args)

    for i in range(1, args.n + 1):
        if time_sampler is not None:
            time_s = sample_timestamp_from_profile(time_rng, time_sampler)
        else:
            time_s = sample_timestamp_offset_seconds(time_rng, args, i)
        if time_s > max_time_s:
            max_time_s = time_s
        frac = time_s / run_duration_s if run_duration_s > 0.0 else 0.0
        read_class = draw_from_mix(rng, class_mix_for_time(frac))
        length = class_length(rng, read_class, args)
        speed_nt_s = class_speed(rng, read_class, args, frac)
        sum_speed += speed_nt_s

        mean_q = class_qscore(rng, read_class, length, speed_nt_s, args, frac)
        barcode = draw_from_mix(rng, args.barcode_mix)
        channel = sample_channel(channel_rng, channel_sampler)
        mux = rng.randint(1, args.muxes)
        is_duplex = 1 if rng.random() < args.duplex_rate else 0
        if is_duplex:
            duplex_n += 1

        class_counts[read_class] = class_counts.get(read_class, 0) + 1
        barcode_counts[barcode] = barcode_counts.get(barcode, 0) + 1
        stats.add(length, mean_q)

        seq = make_seq(rng, length, args.random_seq)
        qchar = phred_char(mean_q)
        qual = qchar * length
        duration_s = float(length) / speed_nt_s

        # Construct one pysam record representation for both SAM and FASTQ output.
        st = iso8601_with_seconds(args.start_dt, time_s)
        duplex_code = 1 if is_duplex else (
            -1 if rng.random() < args.duplex_forming_rate else 0)
        tags = [
            ("RG", read_group_id(args, barcode)),
            ("qs", round(mean_q, 4)),
            ("du", round(duration_s, 6)),
            ("rn", i),
            ("ch", channel),
            ("mx", mux),
            ("st", st),
            ("dx", duplex_code),
        ]
        name = f"{args.id_prefix}_{i:09d}"
        if args.output_format == "sam":
            if rng.random() < args.aligned_rate:
                cigar, nm, reference_span = primary_alignment_cigar(rng, length, mean_q, args)
                reference = "chr2" if rng.random() < args.second_reference_rate else "chr1"
                position = rng.randint(1, args.sam_reference_length - reference_span + 1)
                flag = 16 if rng.random() < args.reverse_rate else 0
                if rng.random() < args.qcfail_rate:
                    flag |= 512
                if rng.random() < args.duplicate_rate:
                    flag |= 1024
                emit_sam_record(
                    sam_writer,
                    sam_records,
                    sam_record(
                        sam_header_value, args, name, flag, reference, position,
                        cigar, seq, qual, [*tags, ("NM", nm)]))
                alignment_counts["primary"] += 1

                if rng.random() < args.secondary_rate:
                    secondary_position = rng.randint(1, args.sam_reference_length - reference_span + 1)
                    emit_sam_record(
                        sam_writer,
                        sam_records,
                        sam_record(
                            sam_header_value, args, name, flag | 256, reference,
                            secondary_position, cigar, seq, qual, [*tags, ("NM", nm)]))
                    alignment_counts["secondary"] += 1
                if rng.random() < args.supplementary_rate:
                    supplementary_position = rng.randint(1, args.sam_reference_length - reference_span + 1)
                    emit_sam_record(
                        sam_writer,
                        sam_records,
                        sam_record(
                            sam_header_value, args, name, flag | 2048, reference,
                            supplementary_position, cigar, seq, qual, [*tags, ("NM", nm)]))
                    alignment_counts["supplementary"] += 1
            else:
                emit_sam_record(
                    sam_writer,
                    sam_records,
                    sam_record(
                        sam_header_value, args, name, 4, "*", 0, None, seq, qual, tags))
                alignment_counts["unmapped"] += 1
        else:
            record = sam_record(
                sam_header_value,
                args,
                name,
                4,
                "*",
                0,
                None,
                seq,
                qual,
                tags)
            write_fastq_record(sys.stdout, record, args, barcode, args.omit_metadata)

        if args.progress_every > 0 and (i % args.progress_every) == 0:
            err.write(f"[simulate-reads] emitted {i} reads\n")

    if sam_records is not None:
        for record in sorted(
                sam_records,
                key=lambda record: (
                    record.reference_id < 0, record.reference_id, record.reference_start)):
            sam_writer.write(record)
    if sam_writer is not None:
        sam_writer.close()

    # Summary to stderr only (safe to pipe FASTQ stdout directly).
    if stats.count > 0:
        mean_len = stats.sum_len / stats.count
        mean_q = stats.sum_q / stats.count
        mean_speed = sum_speed / stats.count
        duration_s = run_duration_s
    else:
        mean_len = 0.0
        mean_q = 0.0
        mean_speed = 0.0
        duration_s = 0.0
        stats.min_len = 0
    err.write(
        "[simulate-reads] "
        f"n={stats.count} seed={args.seed} "
        f"time_profile={args.time_profile} duration_s={duration_s:.1f} max_time_s={max_time_s:.1f} "
        f"len(min/mean/max)={stats.min_len}/{mean_len:.1f}/{stats.max_len} "
        f"q_mean={mean_q:.3f} speed_mean={mean_speed:.2f} "
        f"good={class_counts['good']} mediocre={class_counts['mediocre']} junk={class_counts['junk_short']} duplex={duplex_n}\n"
    )
    if barcode_counts:
        parts = [f"{k}:{v}" for k, v in sorted(barcode_counts.items())]
        err.write(f"[simulate-reads] barcodes {' '.join(parts)}\n")
    if args.output_format == "sam":
        err.write(
            "[simulate-reads] alignments "
            f"primary={alignment_counts['primary']} unmapped={alignment_counts['unmapped']} "
            f"secondary={alignment_counts['secondary']} supplementary={alignment_counts['supplementary']}\n")

    return 0
