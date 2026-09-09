#!/usr/bin/env bash
# Locks in params::validate_alpha_s_running's behavior (see
# docs/PAPER_MAPPING.md, "Known coupling-scheme gaps"): the program must
# reject (alpha_s_running, channel) combinations for which the
# corresponding NLO coefficient function never applies a running-coupling
# factor, rather than silently computing a physically wrong result.
#
# Separate from run.sh's numeric golden-output comparison since these
# cases are about the exit code/rejection itself, not a cross-section
# value -- and they're nearly instant (the guard fires before any physics
# computation), unlike the golden-output cases.
#
# Usage: tests/regression/test_alpha_s_guard.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"
BUILD_DIR="$ROOT_DIR/build"

cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release >/dev/null
cmake --build "$BUILD_DIR" --parallel >/dev/null
BIN="$BUILD_DIR/nlosingleinclusive"

overall_status=0

# name : args : should_reject (1 = must exit non-zero immediately, 0 = must
# NOT exit non-zero within a few seconds -- a genuinely accepted config
# still runs the (slow) physics after the guard passes, so this can't wait
# for it to finish; surviving a short timeout is enough evidence the guard
# didn't reject it)
cases=(
    "smallest+gg:0.05 0.95 0.4 pp 0 g g smallest 10 1:1"
    "smallest+qg:0.05 0.95 0.4 pp 0 q g smallest 10 1:1"
    "smallest+gq:0.05 0.95 0.4 pp 0 g q smallest 10 1:1"
    "daughter+qg:0.05 0.95 0.4 pp 0 q g daughter 10 1:1"
    "daughter+qq:0.05 0.95 0.4 pp 0 q q daughter 10 1:0"
)

for entry in "${cases[@]}"; do
    name="${entry%%:*}"
    rest="${entry#*:}"
    args="${rest%%:*}"
    should_reject="${rest##*:}"

    if [[ "$should_reject" == "1" ]]; then
        if ( cd "$ROOT_DIR" && "$BIN" $args >/tmp/guard_test_out.$$ 2>&1 ); then
            exit_code=0
        else
            exit_code=$?
        fi
        if [[ "$exit_code" -eq 0 ]]; then
            echo "FAIL  $name (expected rejection, program exited 0)" >&2
            overall_status=1
        else
            echo "PASS  $name (rejected, exit $exit_code)"
        fi
    else
        # timeout's own exit code 124 means "still running" -- accepted.
        # Anything else (0 included) after only a few seconds is
        # suspicious for a real physics run, but only a non-zero,
        # non-124 exit indicates the guard rejected it.
        set +e
        ( cd "$ROOT_DIR" && timeout 5 "$BIN" $args >/tmp/guard_test_out.$$ 2>&1 )
        exit_code=$?
        set -e
        if [[ "$exit_code" -ne 0 && "$exit_code" -ne 124 ]]; then
            echo "FAIL  $name (expected acceptance, program exited $exit_code)" >&2
            overall_status=1
        else
            echo "PASS  $name (accepted)"
        fi
    fi
    rm -f /tmp/guard_test_out.$$
done

exit $overall_status
