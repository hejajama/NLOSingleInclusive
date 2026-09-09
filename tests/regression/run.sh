#!/usr/bin/env bash
# Regression harness for the NLOSingleInclusive refactor.
#
# Builds the project, runs each case under tests/regression/cases/<name>/
# (argv in `args`, expected stdout in `golden.csv`), and compares the fresh
# output against the golden file within a numeric tolerance (see
# compare_csv.py). Intended to be run after every stage of the refactor to
# confirm the physics output hasn't changed.
#
# Usage: tests/regression/run.sh [case_name ...]
#   With no arguments, runs every case under cases/.
#   RTOL/ATOL env vars override the default comparison tolerance.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="$ROOT_DIR/build"
CASES_DIR="$SCRIPT_DIR/cases"

RTOL="${RTOL:-1e-6}"
ATOL="${ATOL:-1e-12}"

echo "== configuring/building =="
cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release >/dev/null
cmake --build "$BUILD_DIR" --parallel >/dev/null
BIN="$BUILD_DIR/nlosingleinclusive"

if [[ $# -gt 0 ]]; then
    case_names=("$@")
else
    case_names=()
    for case_path in "$CASES_DIR"/*/; do
        case_names+=("$(basename "$case_path")")
    done
fi

overall_status=0
for name in "${case_names[@]}"; do
    case_dir="$CASES_DIR/$name"
    args_file="$case_dir/args"
    golden_file="$case_dir/golden.csv"

    if [[ ! -f "$args_file" || ! -f "$golden_file" ]]; then
        echo "SKIP  $name (missing args or golden.csv under $case_dir)" >&2
        overall_status=1
        continue
    fi

    args="$(cat "$args_file")"
    actual_file="$(mktemp)"

    # Run from the repo root: the program reads several inputs (BK solution,
    # TA values, ...) via paths relative to the current working directory.
    if ! ( cd "$ROOT_DIR" && "$BIN" $args > "$actual_file" ); then
        echo "FAIL  $name (program exited non-zero)" >&2
        overall_status=1
        rm -f "$actual_file"
        continue
    fi

    if python3 "$SCRIPT_DIR/compare_csv.py" "$golden_file" "$actual_file" \
        --rtol "$RTOL" --atol "$ATOL" --name "$name"; then
        echo "PASS  $name"
    else
        echo "FAIL  $name"
        overall_status=1
    fi
    rm -f "$actual_file"
done

exit $overall_status
