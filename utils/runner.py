#!/usr/bin/env python3
"""Generate a hadron-level p_T spectrum (LO/NLO) for pion production in p+p.

Runs `nlosingleinclusive --level hadron` once per (channel, p_T) pair, saves
each channel's p_T scan to its own file under --tmp-dir, then combines all
channels into a single p_T,LO_total,NLO_total file (--output).

"All channels" means the four (incoming, outgoing) combinations {g,q}x{g,q}
-- g and q (not the individual quark flavors) are enough because
PdfSet::xf/FfSet::zD already sum every light quark+antiquark flavor
internally when incoming/outgoing=="q" (see ../docs/pdf_evaluation_bug.md
for why incoming=="q" didn't always do that correctly).

(channel, p_T) points run concurrently, up to --max-runners at a time; each
individual nlosingleinclusive process still parallelizes its own z-integral
internally over --threads OpenMP threads (see sigma_hadron.cpp), so total
CPU usage is roughly --max-runners x --threads.

Example:
    utils/runner.py --tmp-dir ./runner_tmp --output spectra.csv \\
        --z-points 16 --threads 4 --max-runners 4 \\
        --pt-min 1 --pt-max 10 --pt-step 1 --clean
"""
import argparse
import csv
import os
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

CHANNELS = [("g", "g"), ("g", "q"), ("q", "g"), ("q", "q")]

# Guards stderr progress/warning prints so concurrent runs don't interleave
# mid-line.
_print_lock = threading.Lock()


def log(message):
    with _print_lock:
        print(message, file=sys.stderr)

# Required by the CLI parser even in --level hadron, where they're unused
# (see cli.hpp): zmin is hadron mode's actual z-integral lower bound
# (exposed below as --zmin), zmax/zstep are parton-mode leftovers with no
# effect on hadron-level output.
UNUSED_ZMAX = "0.95"
UNUSED_ZSTEP = "0.5"


def channel_label(incoming, outgoing):
    return f"{incoming}{outgoing}"


def build_command(binary, args, incoming, outgoing, p):
    cmd = [
        binary,
        "--level", "hadron",
        "--col", args.col,
        "--b", str(args.b),
        "--incoming", incoming,
        "--outgoing", outgoing,
        "--rc", args.rc,
        "--p", str(p),
        "--muratio", str(args.muratio),
        "--zmin", str(args.zmin),
        "--zmax", UNUSED_ZMAX,
        "--zstep", UNUSED_ZSTEP,
        "--ff-set", args.ff_set,
        "--z-points", str(args.z_points),
    ]
    if args.sqrts is not None:
        cmd += ["--sqrts", str(args.sqrts)]
    if args.y is not None:
        cmd += ["--y", str(args.y)]
    if args.bk_proton:
        cmd += ["--bk-proton", args.bk_proton]
    if args.bk_nucleus:
        cmd += ["--bk-nucleus", args.bk_nucleus]
    if args.sigma02 is not None:
        cmd += ["--sigma02", str(args.sigma02)]
    return cmd


def run_one(args, incoming, outgoing, p):
    """Runs one (channel, p_T) point. Returns (p_h, LO, NLO), or None (with
    a warning printed) if the run failed or its output couldn't be parsed."""
    cmd = build_command(args.binary, args, incoming, outgoing, p)
    label = channel_label(incoming, outgoing)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(args.threads)

    try:
        result = subprocess.run(cmd, env=env, capture_output=True, text=True,
                                 timeout=args.timeout)
    except subprocess.TimeoutExpired:
        log(f"Warning: channel {label} at p_T={p} timed out after "
            f"{args.timeout}s; skipping this point")
        return None

    if result.returncode != 0:
        log(f"Warning: channel {label} at p_T={p} exited with code "
            f"{result.returncode}; skipping this point. stderr:\n"
            f"{result.stderr.strip()}")
        return None

    lines = [line for line in result.stdout.strip().splitlines() if line.strip()]
    if not lines:
        log(f"Warning: channel {label} at p_T={p} produced no output; "
            "skipping this point")
        return None

    try:
        p_h, lo, nlo = (float(x) for x in lines[-1].split(","))
    except ValueError:
        log(f"Warning: channel {label} at p_T={p} produced unparsable "
            f"output ({lines[-1]!r}); skipping this point")
        return None

    return p_h, lo, nlo


def pt_values(pt_min, pt_max, pt_step):
    values = []
    pt = pt_min
    while pt <= pt_max + 1e-9:
        values.append(round(pt, 10))
        pt += pt_step
    return values


def run_channels(args, pts):
    """Runs every (channel, p_T) point, up to --max-runners at a time, then
    writes each channel's results to its own file under --tmp-dir. Returns
    {label: {p_T: (LO, NLO)}}."""
    os.makedirs(args.tmp_dir, exist_ok=True)
    channel_data = {label: {} for label in
                    (channel_label(i, o) for i, o in CHANNELS)}

    jobs = [(incoming, outgoing, p) for incoming, outgoing in CHANNELS for p in pts]

    with ThreadPoolExecutor(max_workers=args.max_runners) as pool:
        futures = {
            pool.submit(run_one, args, incoming, outgoing, p): (incoming, outgoing, p)
            for incoming, outgoing, p in jobs
        }
        for future in as_completed(futures):
            incoming, outgoing, p = futures[future]
            label = channel_label(incoming, outgoing)
            point = future.result()
            if point is None:
                continue
            p_h, lo, nlo = point
            channel_data[label][p] = (lo, nlo)
            log(f"[{label}] p_T={p} done (LO={lo:.6g}, NLO={nlo:.6g})")

    for incoming, outgoing in CHANNELS:
        label = channel_label(incoming, outgoing)
        channel_file = os.path.join(args.tmp_dir, f"channel_{label}.csv")
        n_ok = 0

        with open(channel_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["p_T", "LO", "NLO"])
            for p in pts:
                point = channel_data[label].get(p)
                if point is None:
                    continue
                lo, nlo = point
                writer.writerow([p, lo, nlo])
                n_ok += 1

        log(f"Channel {label}: {n_ok}/{len(pts)} p_T points succeeded "
            f"-> {channel_file}")
        if n_ok == 0:
            log(f"Warning: channel {label} has no data at all -- it will "
                "be missing from every combined total")

    return channel_data


def combine(args, pts, channel_data):
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["p_T", "LO_total", "NLO_total", "channels_missing"])
        for p in pts:
            lo_sum = 0.0
            nlo_sum = 0.0
            missing = []
            for incoming, outgoing in CHANNELS:
                label = channel_label(incoming, outgoing)
                point = channel_data[label].get(p)
                if point is None:
                    missing.append(label)
                    continue
                lo, nlo = point
                lo_sum += lo
                nlo_sum += nlo

            if missing:
                log(f"Warning: p_T={p} total is missing channel(s) "
                    f"{', '.join(missing)} -- LO_total/NLO_total for this "
                    "point are incomplete")

            writer.writerow([p, lo_sum, nlo_sum, ";".join(missing)])

    log(f"Combined spectrum written to {args.output}")


def clean_tmp(args):
    for incoming, outgoing in CHANNELS:
        channel_file = os.path.join(
            args.tmp_dir, f"channel_{channel_label(incoming, outgoing)}.csv")
        if os.path.exists(channel_file):
            os.remove(channel_file)
    try:
        os.rmdir(args.tmp_dir)
    except OSError:
        pass  # not empty (unrelated files present) -- leave it alone
    log(f"Removed temporary per-channel files under {args.tmp_dir}")


def parse_args():
    ap = argparse.ArgumentParser(
        description="Generate a hadron-level p_T spectrum (LO/NLO) for pion "
                    "production in p+p, summing all four incoming/outgoing "
                    "channels.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    ap.add_argument("--binary", default="./build/nlosingleinclusive",
                     help="path to the built nlosingleinclusive executable")
    ap.add_argument("--tmp-dir", default="./runner_tmp",
                     help="directory the per-channel temporary result files are stored in")
    ap.add_argument("--output", default="spectra.csv",
                     help="file the final combined LO/NLO spectrum is written to")
    ap.add_argument("--z-points", type=int, default=16,
                     help="--z-points forwarded to nlosingleinclusive (Gauss-Legendre "
                          "node count for the hadron-level z integral)")
    ap.add_argument("--threads", type=int, default=1,
                     help="OMP_NUM_THREADS used for each nlosingleinclusive invocation")
    ap.add_argument("--max-runners", type=int, default=1,
                     help="maximum number of (channel, p_T) points run concurrently "
                          "(each still uses --threads OpenMP threads internally, so "
                          "total CPU usage is roughly --max-runners x --threads)")
    ap.add_argument("--clean", action="store_true",
                     help="remove the per-channel temporary files (and --tmp-dir "
                          "itself, if left empty) after combining")

    ap.add_argument("--pt-min", type=float, default=1.0, help="lowest hadron p_T (GeV)")
    ap.add_argument("--pt-max", type=float, default=10.0, help="highest hadron p_T (GeV)")
    ap.add_argument("--pt-step", type=float, default=1.0, help="p_T step (GeV)")

    ap.add_argument("--col", default="pp", help="--col forwarded to nlosingleinclusive")
    ap.add_argument("--b", type=float, default=0.0, help="--b forwarded to nlosingleinclusive")
    ap.add_argument("--rc", default="mom", help="--rc forwarded to nlosingleinclusive")
    ap.add_argument("--muratio", type=float, default=1.0,
                     help="--muratio forwarded to nlosingleinclusive")
    ap.add_argument("--zmin", type=float, default=0.4,
                     help="lower bound of the hadron-level z integral (--zmin)")
    ap.add_argument("--ff-set", default="NNFF10_PIsum_nlo",
                     help="--ff-set forwarded to nlosingleinclusive")
    ap.add_argument("--sqrts", type=float, default=None,
                     help="--sqrts forwarded to nlosingleinclusive (default: the binary's own default)")
    ap.add_argument("--y", type=float, default=None,
                     help="--y forwarded to nlosingleinclusive (default: the binary's own default)")
    ap.add_argument("--bk-proton", default=None, help="--bk-proton forwarded to nlosingleinclusive")
    ap.add_argument("--bk-nucleus", default=None, help="--bk-nucleus forwarded to nlosingleinclusive")
    ap.add_argument("--sigma02", default=None,
                     help="--sigma02 forwarded to nlosingleinclusive (GeV^-2)")
    ap.add_argument("--timeout", type=float, default=1800,
                     help="per-(channel, p_T) run timeout, in seconds")

    return ap.parse_args()


def main():
    args = parse_args()
    pts = pt_values(args.pt_min, args.pt_max, args.pt_step)
    log(f"Running up to {args.max_runners} (channel, p_T) point(s) at a time, "
        f"{args.threads} OpenMP thread(s) each")

    channel_data = run_channels(args, pts)
    combine(args, pts, channel_data)

    if args.clean:
        clean_tmp(args)


if __name__ == "__main__":
    main()
