#!/usr/bin/env python3
"""Generate a hadron-level p_T spectrum (LO/NLO) for pion production in p+p.

Runs `nlosingleinclusive --level hadron` once per (channel, p_T) pair, saves
each channel's p_T scan to its own file under --tmp-dir, then combines all
channels into a single p_T,LO_total,NLO_total,LO_X0_total file (--output).
LO is sigma_LO_k (dipole BK-evolved to Xg, kept for reference/diagnostics),
NLO is the full NLO total (xi-convolution correction plus the X0-frozen LO
baseline), and LO_X0 is that baseline on its own -- the paper's own LO
curve (Eq. 7a/7b evaluated at X0; see ../src/sigma_LO.hpp and
../README.md's "Output"/"Hadron-level cross section" sections).

"All channels" means gg, gq, qg (using the aggregate incoming/outgoing=="q"
token -- PdfSet::xf/FfSet::zD sum every light quark+antiquark flavor
internally for that, see ../docs/pdf_evaluation_bug.md) plus the qq channel,
which -- unlike the other three -- can *not* use incoming=="q" together with
outgoing=="q": that would multiply the flavor-summed incoming PDF by the
flavor-summed outgoing FF, mixing flavors that shouldn't mix (see
../README.md, "Hadron-level flavor sums"; nlosingleinclusive itself refuses
to run that combination). So by default this script instead runs qq once
per explicit matching quark+antiquark flavor (u, d, s, ubar, dbar, sbar) and
sums the six results -- the physically correct sum, at 6x the qq-channel
cost of a single call.

--pi0 trades some of that cost back: if the hadron being fragmented into is
its own antiparticle (like pi0), its FF is the same for a flavor and its
antiflavor (D_u == D_ubar, etc. -- see ../README.md, "pi0 fragmentation
functions"), so f_u*D_u + f_ubar*D_ubar == (f_u+f_ubar)*D_u and each
quark+antiquark pair can be combined into a single nlosingleinclusive call
(--incoming u+ubar --outgoing u) instead of two -- 3 qq-channel calls
instead of 6. This is only correct when --ff-set (or its --pi0 default) is
actually charge-conjugation symmetric like that; it is not a safe default
for an arbitrary (non-self-conjugate) hadron.

(channel, p_T) points run concurrently, up to --max-runners at a time; each
individual nlosingleinclusive process still parallelizes its own z-integral
internally over --threads OpenMP threads (see sigma_hadron.cpp), so total
CPU usage is roughly --max-runners x --threads.

Example:
    utils/runner.py --tmp-dir ./runner_tmp --output spectra.csv \\
        --z-points 16 --threads 4 --max-runners 4 \\
        --pt-min 1 --pt-max 10 --pt-step 1 --clean

    # pi0 production, using the faster qq pairing:
    utils/runner.py --pi0 --tmp-dir ./runner_tmp --output spectra_pi0.csv \\
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

BASE_CHANNELS = [("g", "g"), ("g", "q"), ("q", "g")]

QQ_FLAVORS = ["u", "d", "s", "ubar", "dbar", "sbar"]

# Default qq-channel treatment: one call per explicit matching flavor,
# summed -- correct for any hadron, no assumption about its FF.
QQ_CHANNELS_FULL = [(f, f) for f in QQ_FLAVORS]

# --pi0's qq-channel treatment: one call per quark+antiquark pair, using
# PdfSet::xf's '+'-joined incoming-flavor sum (see pdf_set.cpp) -- correct
# only because D_u^pi0==D_ubar^pi0 (and same for d/s), so pairing them
# behind a single explicit outgoing flavor ("u" stands in for both u and
# ubar's FF) reproduces the sum of running them separately, in one
# nlosingleinclusive call's worth of time instead of two -- see the module
# docstring above and ../README.md, "pi0 fragmentation functions".
QQ_CHANNELS_PI0 = [("u+ubar", "u"), ("d+dbar", "d"), ("s+sbar", "s")]

# --ff-set default used when --pi0 is given and --ff-set is not: the
# isospin/charge-conjugation-symmetric pi0 FF this script (and
# nlosingleinclusive's own --ff-set comma convention, see ../README.md)
# builds as the average of the charged-pion sets.
DEFAULT_PI0_FF_SET = "NNFF10_PIp_nlo,NNFF10_PIm_nlo"


def channels_for(args):
    return BASE_CHANNELS + (QQ_CHANNELS_PI0 if args.pi0 else QQ_CHANNELS_FULL)

# Guards stderr progress/warning prints so concurrent runs don't interleave
# mid-line.
_print_lock = threading.Lock()


def log(message):
    with _print_lock:
        print(message, file=sys.stderr)


def channel_label(incoming, outgoing):
    return f"{incoming}{outgoing}"


def build_command(binary, args, incoming, outgoing, p):
    cmd = [
        binary,
        "--level", "hadron",
        "--col", args.col,
        "--incoming", incoming,
        "--outgoing", outgoing,
        "--rc", args.rc,
        "--pt", str(p),
        "--muratio", str(args.muratio),
        "--zmin", str(args.zmin),
        "--z-points", str(args.z_points),
    ]
    if args.b is not None:
        cmd += ["--b", str(args.b)]
    if args.ff_set is not None:
        cmd += ["--ff-set", args.ff_set]
    if args.pdf_set is not None:
        cmd += ["--pdf-set", args.pdf_set]
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
    """Runs one (channel, p_T) point. Returns (p_h, LO, NLO, LO_X0), or None
    (with a warning printed) if the run failed or its output couldn't be
    parsed. LO is sigma_LO_k (dipole BK-evolved to Xg, kept for reference),
    NLO is the full NLO total (xi-convolution correction + the X0-frozen LO
    baseline), and LO_X0 is that X0 baseline on its own -- the paper's own
    LO curve (Eq. 7a/7b at X0; see src/sigma_LO.hpp's sigma_LO_k_X0 doc
    comment and README.md's "Output"/"Hadron-level cross section" sections)."""
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
        p_h, lo, nlo, lo_x0 = (float(x) for x in lines[-1].split(","))
    except ValueError:
        log(f"Warning: channel {label} at p_T={p} produced unparsable "
            f"output ({lines[-1]!r}); skipping this point")
        return None

    return p_h, lo, nlo, lo_x0


def pt_values(pt_min, pt_max, pt_step):
    values = []
    pt = pt_min
    while pt <= pt_max + 1e-9:
        values.append(round(pt, 10))
        pt += pt_step
    return values


def run_channels(args, pts, channels):
    """Runs every (channel, p_T) point, up to --max-runners at a time, then
    writes each channel's results to its own file under --tmp-dir. Returns
    {label: {p_T: (LO, NLO, LO_X0)}}."""
    os.makedirs(args.tmp_dir, exist_ok=True)
    channel_data = {label: {} for label in
                    (channel_label(i, o) for i, o in channels)}

    jobs = [(incoming, outgoing, p) for incoming, outgoing in channels for p in pts]

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
            p_h, lo, nlo, lo_x0 = point
            channel_data[label][p] = (lo, nlo, lo_x0)
            log(f"[{label}] p_T={p} done (LO={lo:.6g}, NLO={nlo:.6g}, LO_X0={lo_x0:.6g})")

    for incoming, outgoing in channels:
        label = channel_label(incoming, outgoing)
        channel_file = os.path.join(args.tmp_dir, f"channel_{label}.csv")
        n_ok = 0

        with open(channel_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["p_T", "LO", "NLO", "LO_X0"])
            for p in pts:
                point = channel_data[label].get(p)
                if point is None:
                    continue
                lo, nlo, lo_x0 = point
                writer.writerow([p, lo, nlo, lo_x0])
                n_ok += 1

        log(f"Channel {label}: {n_ok}/{len(pts)} p_T points succeeded "
            f"-> {channel_file}")
        if n_ok == 0:
            log(f"Warning: channel {label} has no data at all -- it will "
                "be missing from every combined total")

    return channel_data


def combine(args, pts, channel_data, channels):
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["p_T", "LO_total", "NLO_total", "LO_X0_total", "channels_missing"])
        for p in pts:
            lo_sum = 0.0
            nlo_sum = 0.0
            lo_x0_sum = 0.0
            missing = []
            for incoming, outgoing in channels:
                label = channel_label(incoming, outgoing)
                point = channel_data[label].get(p)
                if point is None:
                    missing.append(label)
                    continue
                lo, nlo, lo_x0 = point
                lo_sum += lo
                nlo_sum += nlo
                lo_x0_sum += lo_x0

            if missing:
                log(f"Warning: p_T={p} total is missing channel(s) "
                    f"{', '.join(missing)} -- LO_total/NLO_total/LO_X0_total "
                    "for this point are incomplete")

            writer.writerow([p, lo_sum, nlo_sum, lo_x0_sum, ";".join(missing)])

    log(f"Combined spectrum written to {args.output}")


def clean_tmp(args, channels):
    for incoming, outgoing in channels:
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
                    "production in p+p, summing the gg/gq/qg channels and "
                    "every qq-channel quark flavor (or, with --pi0, each "
                    "quark+antiquark pair).",
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
    ap.add_argument("--pi0", action="store_true",
                     help="fast path for pi0 production: runs the qq channel as 3 "
                          "quark+antiquark-pair calls (u+ubar, d+dbar, s+sbar) instead of "
                          "6 single-flavor calls, exploiting D_u^pi0==D_ubar^pi0 (and d/s) "
                          "-- roughly halves qq-channel runtime. Only correct when --ff-set "
                          "is charge-conjugation symmetric like a pi0 FF; if --ff-set isn't "
                          f"given, defaults it to '{DEFAULT_PI0_FF_SET}' instead of the "
                          "binary's own default. See ../README.md, \"pi0 fragmentation "
                          "functions\"")

    ap.add_argument("--pt-min", type=float, default=1.0, help="lowest hadron p_T (GeV)")
    ap.add_argument("--pt-max", type=float, default=10.0, help="highest hadron p_T (GeV)")
    ap.add_argument("--pt-step", type=float, default=1.0, help="p_T step (GeV)")

    ap.add_argument("--col", default="pp", help="--col forwarded to nlosingleinclusive")
    ap.add_argument("--b", type=float, default=None,
                     help="--b forwarded to nlosingleinclusive (impact parameter; only "
                          "meaningful, and required by the binary, for --col pA -- omit "
                          "it for the default --col pp)")
    ap.add_argument("--rc", default="mom", help="--rc forwarded to nlosingleinclusive")
    ap.add_argument("--muratio", type=float, default=1.0,
                     help="--muratio forwarded to nlosingleinclusive")
    ap.add_argument("--zmin", type=float, default=0.4,
                     help="lower bound of the hadron-level z integral (--zmin)")
    ap.add_argument("--ff-set", default=None,
                     help="--ff-set forwarded to nlosingleinclusive (default: the binary's "
                          f"own default, unless --pi0 is given, in which case '{DEFAULT_PI0_FF_SET}')")
    ap.add_argument("--pdf-set", default=None,
                     help="--pdf-set forwarded to nlosingleinclusive (default: the binary's own default, params::pdfname)")
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
    if args.pi0 and args.ff_set is None:
        args.ff_set = DEFAULT_PI0_FF_SET
        log(f"--pi0 given without --ff-set: defaulting --ff-set to {DEFAULT_PI0_FF_SET}")

    channels = channels_for(args)
    pts = pt_values(args.pt_min, args.pt_max, args.pt_step)
    log(f"Running up to {args.max_runners} (channel, p_T) point(s) at a time, "
        f"{args.threads} OpenMP thread(s) each")

    channel_data = run_channels(args, pts, channels)
    combine(args, pts, channel_data, channels)

    if args.clean:
        clean_tmp(args, channels)


if __name__ == "__main__":
    main()
