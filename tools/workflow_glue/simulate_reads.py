"""Stream synthetic FASTQ records to stdout.

Designed for fastcat testing:
- read lengths sampled from a Gamma distribution
- bimodal mean-q distribution (two Gaussian peaks, e.g. ~10 and ~20)
- explicit SAM-style `qs:f:` tag in FASTQ comment so fastcat can use it
  directly without recomputing quality

Example:
  python3 tools/dev-tools simulate_reads \
    --n 200000 --shape 2.6 --scale 30000 \
    --q1 10 --q2 20 --q1-weight 0.6 \
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


DEFAULT_RUNID = "12345678-1234-1234-1234-123456789abc"
DEFAULT_MODEL = "dna_r10.4.1_e8.2_400bps_hac@v5.0.0"
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


def sample_qscore(
    rng,
    q1,
    q2,
    q1_sd,
    q2_sd,
    q1_weight,
    qmin,
    qmax,
):
    peak = 1 if rng.random() < q1_weight else 2
    if peak == 1:
        q = rng.gauss(q1, q1_sd)
    else:
        q = rng.gauss(q2, q2_sd)
    q = clamp(q, qmin, qmax)
    return q, peak


def sample_length(
    rng,
    shape,
    scale,
    min_len,
    max_len,
):
    length = int(round(rng.gammavariate(shape, scale)))
    if length < min_len:
        return min_len
    if length > max_len:
        return max_len
    return length


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


def argparser():
    p = argparse.ArgumentParser(
        description="Stream synthetic FASTQ with gamma lengths and bimodal mean-q; emits qs:f tag.",
        add_help=False,
    )

    run = p.add_argument_group("simulation")
    run.add_argument("--n", type=int, default=10000, help="Number of reads to emit (default: 10000)")
    run.add_argument("--seed", type=int, default=42, help="RNG seed (default: 42)")

    length = p.add_argument_group("read length distribution")
    length.add_argument("--shape", type=float, default=2.6, help="Gamma shape k (default: 2.6)")
    length.add_argument("--scale", type=float, default=30000.0, help="Gamma scale theta (default: 30000)")
    length.add_argument("--min-len", type=int, default=100, help="Minimum read length clamp (default: 100)")
    length.add_argument("--max-len", type=int, default=1_000_000, help="Maximum read length clamp (default: 1000000)")

    qscore = p.add_argument_group("qscore distribution")
    qscore.add_argument("--q1", type=float, default=10.0, help="First qscore peak mean (default: 10)")
    qscore.add_argument("--q2", type=float, default=20.0, help="Second qscore peak mean (default: 20)")
    qscore.add_argument("--q1-sd", type=float, default=1.6, help="Stddev around q1 peak (default: 1.6)")
    qscore.add_argument("--q2-sd", type=float, default=1.8, help="Stddev around q2 peak (default: 1.8)")
    qscore.add_argument("--q1-weight", type=float, default=0.6, help="Mix weight for q1 peak in [0,1] (default: 0.6)")
    qscore.add_argument("--q-min", type=float, default=1.0, help="Minimum qscore clamp (default: 1)")
    qscore.add_argument("--q-max", type=float, default=40.0, help="Maximum qscore clamp (default: 40)")

    barcode = p.add_argument_group("barcode mixture")
    barcode.add_argument("--barcodes", nargs="+", default=["unclassified"], help="Barcode names to sample")
    barcode.add_argument("--barcode-weights", nargs="+", type=float, default=[1.0], help="Weights matching --barcodes")

    metadata = p.add_argument_group("read metadata")
    metadata.add_argument("--runid", default=DEFAULT_RUNID, help=f"Run ID for RG composition (default: {DEFAULT_RUNID})")
    metadata.add_argument("--model", default=DEFAULT_MODEL, help=f"Basecaller model in RG (default: {DEFAULT_MODEL})")
    metadata.add_argument("--id-prefix", default="simread", help="Read ID prefix (default: simread)")
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
        help="Timestamp spacing profile (default: saturating)")
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

    output = p.add_argument_group("output")
    output.add_argument("--random-seq", action="store_true", help="Emit random A/C/G/T instead of homopolymer A")
    output.add_argument("--progress-every", type=int, default=0, help="If >0, print progress to stderr every N reads")

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
    if args.progress_every < 0:
        p.error("--progress-every must be >= 0")

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


def main(args):
    validate_args(argparser(), args)
    rng = random.Random(args.seed)
    time_rng = random.Random(args.seed + 1_000_003)
    channel_rng = random.Random(args.seed + 2_000_003)
    channel_sampler = build_channel_sampler(args.device, args.channels)

    out = sys.stdout
    err = sys.stderr

    stats = RunningStats()
    qpeak1 = 0
    qpeak2 = 0
    duplex_n = 0
    max_time_s = 0.0
    barcode_counts = {name: 0 for name, _ in args.barcode_mix}

    for i in range(1, args.n + 1):
        length = sample_length(rng, args.shape, args.scale, args.min_len, args.max_len)
        mean_q, peak = sample_qscore(
            rng,
            args.q1,
            args.q2,
            args.q1_sd,
            args.q2_sd,
            args.q1_weight,
            args.q_min,
            args.q_max,
        )
        barcode = draw_from_mix(rng, args.barcode_mix)
        channel = sample_channel(channel_rng, channel_sampler)
        mux = rng.randint(1, args.muxes)
        is_duplex = 1 if rng.random() < args.duplex_rate else 0
        if is_duplex:
            duplex_n += 1

        if peak == 1:
            qpeak1 += 1
        else:
            qpeak2 += 1
        barcode_counts[barcode] = barcode_counts.get(barcode, 0) + 1
        stats.add(length, mean_q)

        seq = make_seq(rng, length, args.random_seq)
        qchar = phred_char(mean_q)
        qual = qchar * length

        # Keep tags SAM-like so fastcat ingest reads them as BAM aux fields directly.
        time_s = sample_timestamp_offset_seconds(time_rng, args, i)
        if time_s > max_time_s:
            max_time_s = time_s
        st = iso8601_with_seconds(args.start_dt, time_s)
        rg = f"{args.runid}_{args.model}"
        if barcode and barcode != "unclassified":
            rg = f"{rg}_{barcode}"
        tags = (
            f"qs:f:{mean_q:.4f}\t"
            f"rn:i:{i}\t"
            f"ch:i:{channel}\t"
            f"mx:i:{mux}\t"
            f"st:Z:{st}\t"
            f"SM:Z:{barcode}\t"
            f"RG:Z:{rg}\t"
            f"dx:i:{is_duplex}"
        )

        out.write(f"@{args.id_prefix}_{i:09d} {tags}\n")
        out.write(f"{seq}\n+\n{qual}\n")

        if args.progress_every > 0 and (i % args.progress_every) == 0:
            err.write(f"[simulate-reads] emitted {i} reads\n")

    # Summary to stderr only (safe to pipe FASTQ stdout directly).
    if stats.count > 0:
        mean_len = stats.sum_len / stats.count
        mean_q = stats.sum_q / stats.count
        duration_s = timestamp_duration_seconds(args)
    else:
        mean_len = 0.0
        mean_q = 0.0
        duration_s = 0.0
        stats.min_len = 0
    err.write(
        "[simulate-reads] "
        f"n={stats.count} seed={args.seed} "
        f"time_profile={args.time_profile} duration_s={duration_s:.1f} max_time_s={max_time_s:.1f} "
        f"len(min/mean/max)={stats.min_len}/{mean_len:.1f}/{stats.max_len} "
        f"q_mean={mean_q:.3f} qpeak1={qpeak1} qpeak2={qpeak2} duplex={duplex_n}\n"
    )
    if barcode_counts:
        parts = [f"{k}:{v}" for k, v in sorted(barcode_counts.items())]
        err.write(f"[simulate-reads] barcodes {' '.join(parts)}\n")

    return 0
