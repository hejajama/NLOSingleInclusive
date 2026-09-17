#!/usr/bin/env python3
"""Tolerance-based CSV comparison for the NLOSingleInclusive regression suite.

The program under test prints floating-point results with only ~6 significant
digits, so a byte-exact diff is too strict to survive things like refactoring
the order of arithmetic operations that are mathematically but not bit-for-bit
equivalent. Instead this compares row-by-row, column-by-column, with the
classic isclose tolerance:

    abs(actual - golden) <= atol + rtol * abs(golden)

Exit code 0 = match within tolerance, 1 = mismatch (row/column count or
value), 2 = usage/IO error.
"""
import argparse
import sys


def load_csv(path):
    rows = []
    with open(path, newline="") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                rows.append([float(x) for x in line.split(",")])
            except ValueError:
                continue  # header row (e.g. "z,k [GeV],...")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("golden", help="path to the golden/reference CSV")
    ap.add_argument("actual", help="path to the freshly produced CSV")
    ap.add_argument("--rtol", type=float, default=1e-6, help="relative tolerance")
    ap.add_argument("--atol", type=float, default=1e-12, help="absolute tolerance floor")
    ap.add_argument("--name", default="", help="case name, used only in messages")
    args = ap.parse_args()

    label = f"[{args.name}] " if args.name else ""

    try:
        golden = load_csv(args.golden)
        actual = load_csv(args.actual)
    except (OSError, ValueError) as exc:
        print(f"{label}error reading input: {exc}", file=sys.stderr)
        return 2

    if len(golden) != len(actual):
        print(
            f"{label}row count mismatch: golden={len(golden)} actual={len(actual)}",
            file=sys.stderr,
        )
        return 1

    failures = []
    max_rel = 0.0
    for i, (grow, arow) in enumerate(zip(golden, actual)):
        if len(grow) != len(arow):
            failures.append(
                f"row {i}: column count mismatch (golden={len(grow)}, actual={len(arow)})"
            )
            continue
        for j, (g, a) in enumerate(zip(grow, arow)):
            tol = args.atol + args.rtol * abs(g)
            diff = abs(a - g)
            if diff > tol:
                failures.append(
                    f"row {i} col {j}: golden={g!r} actual={a!r} diff={diff:.3e} > tol={tol:.3e}"
                )
            if abs(g) > 0:
                max_rel = max(max_rel, diff / abs(g))

    if failures:
        print(f"{label}{len(failures)} mismatch(es) (showing up to 10):", file=sys.stderr)
        for line in failures[:10]:
            print(f"  {line}", file=sys.stderr)
        return 1

    print(f"{label}OK ({len(golden)} rows, max relative diff {max_rel:.3e})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
