# NLOSingleInclusive

Parton-level LO + NLO single-inclusive cross section for forward
particle production in proton-proton and proton-nucleus collisions in
the Color Glass Condensate framework.

## Reference

> Heikki Mäntysaari and Yossathorn Tawabutr, *"Complete Next-to-Leading
> Order Calculation of Single Inclusive π⁰ Production in Forward
> Proton-Nucleus Collisions"*, [arXiv:2310.06640](https://arxiv.org/abs/2310.06640).

See [docs/PAPER_MAPPING.md](docs/PAPER_MAPPING.md) for how the paper's
equations map onto the code.

## What it computes

For a fixed transverse momentum `p` (the produced parton's `p_T`) and
fixed rapidity `y`, the program scans `z = p/k` (`k` being the incoming
parton's transverse momentum, so `k = p/z`) and prints the LO and NLO
partonic cross sections at each point (Sec. IV of the paper). This is
the parton-level result only -- no fragmentation-function convolution to
the hadron level (Sec. V) is performed.

The incoming/outgoing parton channel (qq, qg, gq or gg), the
running-coupling prescription, the collision system (`pp` or `pA`), and
the impact parameter `b` (for `pA`) are all chosen on the command line.
Everything else (BK-solution file, PDF set, `sigma0`, `Qs0^2`, ...) is a
compile-time constant in [src/params.hpp](src/params.hpp).

## Building

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

Requires GSL and LHAPDF (see [CMakeLists.txt](CMakeLists.txt); set
`-DLHAPDF_ROOT=/path/to/lhapdf` if it isn't auto-detected).

## Running

```sh
./build/nlosingleinclusive \
    --zmin 0.05 --zmax 0.95 --zstep 0.4 \
    --col pp --b 0 \
    --incoming q --outgoing q \
    --rc mom \
    --p 10 --muratio 1
```

Flags (any order, all required unless marked optional):

| Flag | Meaning |
|---|---|
| `--zmin`, `--zmax`, `--zstep` | Scan range/step in `z = p/k` |
| `--col` | Collision system: `pp` or `pA` |
| `--b` | Impact parameter (fm); looked up in the `TAvalues` file for `pA` |
| `--incoming` | Incoming parton: `g`, `q`, `u`, `d`, `s`, `ubar`, `dbar`, `sbar` |
| `--outgoing` | Outgoing parton: `q` or `g` (selects the qq/qg/gq/gg channel together with `--incoming`) |
| `--rc` | Running-coupling prescription: `fixed`, `mom`, `parent`, `daughter`, `smallest`, `mixed`, `mixedbd` |
| `--p` | Produced parton's transverse momentum (GeV) |
| `--muratio` | `mu / p` ratio setting the factorization/renormalization scale, `mu^2 = (muratio * p)^2` |
| `--sqrts` *(optional)* | Collision energy sqrt(s) (GeV); defaults to `params::SQRTS` in [src/params.hpp](src/params.hpp) |
| `--y` *(optional)* | Produced parton's rapidity; defaults to `params::yh` in [src/params.hpp](src/params.hpp) |
| `--bk-proton` *(optional)* | BK solution file used when `--col pp`; defaults to `params::bksolpp` |
| `--bk-nucleus` *(optional)* | BK solution filename *prefix* used when `--col pA`; defaults to `params::bksolpA`. The code appends the digits of `--b` directly (no separator inserted), so a prefix meant to read as `..._<b>` must already end in `_`, e.g. `--bk-nucleus ./KCBK_fit_3/Pb_b_` |

A legacy positional form is also still accepted for backward
compatibility with older scripts:

```sh
./build/nlosingleinclusive zmin zmax zstep col b incoming outgoing rc p muratio
```

(same fields, in that fixed order; `--sqrts`/`--y`/`--bk-proton`/
`--bk-nucleus` have no positional equivalent, so that form always uses the
compiled-in defaults).

## Output

One CSV line per scanned `z`, with no header, in decreasing `z` from
`zmax` to `zmin`:

```
z,k,sigma_LO,sigma_NLO
```

- `z` -- the scan variable, `z = p/k`
- `k` -- the incoming parton's transverse momentum (GeV), `k = p/z`
- `sigma_LO` -- the leading-order partonic cross section at this point
- `sigma_NLO` -- the next-to-leading-order correction at this point

e.g.:

```
0.95,10.5263,-2.94206e-07,4.43243e-07
0.55,18.1818,-3.32203e-08,1.27906e-08
0.15,66.6667,-7.92952e-12,-3.70907e-11
```
