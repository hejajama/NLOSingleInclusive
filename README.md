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
| `--b` | Impact parameter ($\mathrm{GeV}^{-1}$); looked up in the `TAvalues` file for `pA` |
| `--incoming` | Incoming parton: `g`, `q`, `u`, `d`, `s`, `ubar`, `dbar`, `sbar` |
| `--outgoing` | Outgoing parton: `g`, `q`, `u`, `d`, `s`, `ubar`, `dbar`, `sbar` (selects the qq/qg/gq/gg channel together with `--incoming`; `g`/non-`g` is what determines the channel, so e.g. `--outgoing u` selects the same channel as `--outgoing q`) |
| `--rc` | Running-coupling prescription: `fixed`, `mom`, `parent`, `daughter`, `smallest`, `mixed`, `mixedbd` |
| `--p` | Produced parton's transverse momentum (GeV) |
| `--muratio` | `mu / p` ratio setting the factorization/renormalization scale, `mu^2 = (muratio * p)^2` |
| `--sqrts` *(optional)* | Collision energy sqrt(s) (GeV); defaults to `params::SQRTS` in [src/params.hpp](src/params.hpp) |
| `--y` *(optional)* | Produced parton's rapidity; defaults to `params::yh` in [src/params.hpp](src/params.hpp) |
| `--bk-proton` *(optional)* | BK solution file used when `--col pp`; defaults to `params::bksolpp` |
| `--bk-nucleus` *(optional)* | BK solution filename *prefix* used when `--col pA`; defaults to `params::bksolpA`. The code appends the digits of `--b` directly (no separator inserted), so a prefix meant to read as `..._<b>` must already end in `_`, e.g. `--bk-nucleus ./KCBK_fit_3/Pb_b_` |
| `--sigma02` *(optional)* | `sigma0/2`, in GeV^-2, doubled into `sigma0` (Sr_0's pA-branch normalization, Eq. 14); defaults to `params::sigma0/2`. Append a literal `mb` to give the value in millibarns instead, e.g. `--sigma02 47.2 mb` |

A legacy positional form is also still accepted for backward
compatibility with older scripts:

```sh
./build/nlosingleinclusive zmin zmax zstep col b incoming outgoing rc p muratio
```

(same fields, in that fixed order; `--sqrts`/`--y`/`--bk-proton`/
`--bk-nucleus`/`--sigma02` have no positional equivalent, so that form
always uses the compiled-in defaults).

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

## Hadron-level cross section (`--level hadron`)

Passing `--level hadron` instead convolves the same parton-level LO/NLO
cross section with a collinear fragmentation function `D(z,Q^2)` over
`z = p_h/k` (Sec. V), reporting a single result at the fixed hadron
transverse momentum `--pt` = `p_h` (see
[src/sigma_hadron.hpp](src/sigma_hadron.hpp)):

```sh
./build/nlosingleinclusive \
    --level hadron --zmin 0.4 \
    --col pp --b 0 \
    --incoming u --outgoing u \
    --rc mom --pt 2 --muratio 1
```

`--zmax`/`--zstep` are not needed for `--level hadron` (its z integral
only uses `--zmin`; the upper bound is always 1). A few extra flags apply
only to this mode:

| Flag | Meaning |
|---|---|
| `--ff-set` *(optional)* | LHAPDF fragmentation-function set (default `params::default_ffname`). Can also be two set names separated by a comma (no spaces) to average them -- see "pi0 fragmentation functions" below |
| `--pdf-set` *(optional)* | LHAPDF PDF set for the incoming parton (default `params::pdfname`); also usable with `--level parton` |
| `--z-points` *(optional)* | Gauss-Legendre node count for the z integral (default 16; see the validation note in [src/sigma_hadron.cpp](src/sigma_hadron.cpp)) |

### Output

A single CSV line (plus header):

```
pT [GeV],dN_LO(x_g) / d^2p_Tdy [1/GeV^2],dN_NLO / d^2p_Tdy [1/GeV^2],dN_LO(X0) / d^2p_Tdy [1/GeV^2]
```

(the header instead reads `d sigma_LO(x_g) / ...`, `d sigma_NLO / ...`,
`d sigma_LO(X0) / ...` in `[1/GeV^4]` for `--col pp`, where the output is
already multiplied by `sigma0/2`.) The three value columns are:

- `LO(x_g)` -- `sigma_LO_k` (Eq. 7a/7b), with the dipole amplitude
  BK-evolved to `Xg` (Eq. 5). This does **not** match how the paper itself
  defines the LO curve (see `LO(X0)` below) -- kept only for
  reference/diagnostics.
- `NLO` -- the full NLO-accuracy result: the `xi`-convolution correction
  (Eqs. 9/11a-11c) plus `LO(X0)`. This is the paper's own NLO curve.
- `LO(X0)` -- Eq. 7a/7b evaluated at the dipole's initial-condition scale
  `X0` (Eq. 13/14), exactly as the paper's text specifies right before
  Eq. (7a)/(7b): *"dipole amplitudes in the LO term are evaluated at the
  initial scale, X0, of the small-x evolution"* -- the rapidity evolution
  down to `Xg` only enters through the NLO `xi`-convolution's `X(xi)`
  (Eq. 4). This is the paper's own LO curve.

See [src/sigma_hadron.hpp](src/sigma_hadron.hpp)'s `HadronSigma` and
[src/sigma_LO.hpp](src/sigma_LO.hpp)'s `sigma_LO_k`/`sigma_LO_k_X0` for
the corresponding code.

### Hadron-level flavor sums

`--incoming q` and `--outgoing q` each sum the corresponding PDF/FF over
all `Nf` light quark *and* antiquark flavors (`PdfSet::xf`/`FfSet::zD`).
At parton level (`--level parton`, no fragmentation function involved)
that's fine on its own: each incoming flavor contributes to the
inclusive quark yield independently, so summing the PDF alone is
correct.

At hadron level, **a multi-flavor `--incoming` (`q`, or a `+`-joined token
like `u+ubar`, see below) together with `--outgoing q` is wrong**, and
the program refuses to run it (see the check at the top of
`sigma_hadron_ph` in [src/sigma_hadron.cpp](src/sigma_hadron.cpp)). The
qq channel's hard function conserves quark flavor -- an incoming quark
of flavor `i` stays flavor `i` all the way through to the fragmenting
parton -- so the physically correct hadron-level sum is
`sum_i f_i(x) * D_i(z)`, one term per matching flavor pair. Aggregating
both sides instead gives `(sum_i f_i) * (sum_j D_j)`, which includes
unphysical cross-flavor terms (e.g. an incoming u-quark weighted by the
d-quark FF).

**Preferred approach:** to get the full flavor-summed qq hadron-level
cross section, run once per explicit matching flavor pair and add the
resulting LO/NLO values:

```
--incoming u    --outgoing u
--incoming d    --outgoing d
--incoming s    --outgoing s
--incoming ubar --outgoing ubar
--incoming dbar --outgoing dbar
--incoming sbar --outgoing sbar
```

(`--outgoing` accepts these same explicit flavor tokens as `--incoming`,
mirroring `PdfSet::xf`, precisely so this per-flavor loop is possible.)

`--incoming q --outgoing g` (qg) and `--incoming g --outgoing q` (gq) do
*not* have this problem and can keep using the aggregate `q` token as-is:
only one side (the PDF for qg, the FF for gq) ever sums over flavor in
those channels, multiplied by a flavor-independent gluon PDF/FF on the
other side, so summing first and multiplying after gives the same result
as summing the per-flavor products. `--incoming g --outgoing g` (gg) has
no quark-flavor sum on either side at all. `utils/runner.py` (the
`p_T`-spectrum driver) implements exactly this: gg/gq/qg each run once
with the aggregate token, while qq runs once per explicit matching
flavor pair (see below) and sums the six results.

`--incoming` also accepts several flavor tokens joined by `+` (e.g.
`u+ubar`), summing just those flavors' PDFs rather than all `Nf` of
them (`PdfSet::xf`, [src/pdf_set.cpp](src/pdf_set.cpp)) -- this is safe
to combine with a *single* explicit `--outgoing` flavor (it's only
`--outgoing q`, the all-flavor FF sum, that the check above rejects),
and is how the `--pi0` fast path below halves the qq-channel cost.

### pi0 fragmentation functions

LHAPDF ships no charge-neutral pion (pi0) fragmentation-function set --
only charged-pion sets such as `NNFF10_PIp_nlo` (pi+), `NNFF10_PIm_nlo`
(pi-), and `NNFF10_PIsum_nlo` (the pi+ + pi- *charged*-pion sum, not
pi0). The standard workaround, used here, is isospin symmetry:
`D_pi0(z,Q^2) = (D_pi+(z,Q^2) + D_pi-(z,Q^2)) / 2`.

`--ff-set` supports this directly: pass two LHAPDF set names separated by
a comma (no spaces) and `FfSet::zD` averages `z*D(z,Q^2)` over both sets
(see [src/ff_set.hpp](src/ff_set.hpp)/[.cpp](src/ff_set.cpp)) --

```sh
./build/nlosingleinclusive \
    --level hadron --zmin 0.4 \
    --col pp --b 0 \
    --incoming u --outgoing u \
    --rc mom --pt 2 --muratio 1 \
    --ff-set NNFF10_PIp_nlo,NNFF10_PIm_nlo
```

`FfSet::name()` returns the comma-joined string unchanged, so this also
works transparently with `sigma_hadron.cpp`'s per-OpenMP-thread `FfSet`
copies (each thread rebuilds its own `FfSet` from `ff.name()`) -- no
special-casing needed there. The two set names can be anything LHAPDF
recognizes; this isn't pi0-specific, just the convention used to get a
pi0 FF out of two charged-pion sets.

### `utils/runner.py --pi0`: a faster qq channel for pi0

pi0 is its own antiparticle, so its FF is the same for a quark and its
antiquark: `D_u^pi0(z) = D_ubar^pi0(z)` (and likewise for d/s), by
charge-conjugation symmetry combined with the isospin averaging above --
crossing quark <-> antiquark and hadron <-> antihadron leaves a
self-conjugate hadron's FF unchanged. That symmetry lets the qq-channel
sum for a matching quark+antiquark pair be computed as a *single*
`nlosingleinclusive` call instead of two:

```
f_u(x)*D_u(z) + f_ubar(x)*D_ubar(z) = [f_u(x) + f_ubar(x)] * D_u(z)
```

This is exact, not an approximation, given `sigma_LO_k`/`sigma_NLO_k` are
linear in the incoming PDF for the QQ channel (`PointTables`'s expensive
part -- `build_coefficient_tables`'s NLO-coefficient double integrals,
[src/point_tables.cpp](src/point_tables.cpp) -- doesn't depend on which
quark flavor is incoming at all, only on the channel; flavor only enters
through the cheap PDF lookups in `build_xi_convolution`). So a single
`--incoming u+ubar --outgoing u` call (using `PdfSet::xf`'s `+`-joined
incoming sum above) reproduces the sum of separate `--incoming u` and
`--incoming ubar` runs, while building `PointTables` -- the dominant cost
-- only once instead of twice.

`utils/runner.py --pi0` uses exactly this: the qq channel runs as 3 calls
(`u+ubar`/`u`, `d+dbar`/`d`, `s+sbar`/`s`) instead of 6, roughly halving
qq-channel runtime, and defaults `--ff-set` to
`NNFF10_PIp_nlo,NNFF10_PIm_nlo` if `--ff-set` isn't given explicitly.
This is only correct because the configured FF really is
charge-conjugation symmetric -- `--pi0` (and the `u+ubar`-style incoming
token generally) is not a safe default for an arbitrary, non-self-
conjugate hadron, where a flavor and its antiflavor's FFs genuinely
differ.

```sh
utils/runner.py --pi0 --tmp-dir ./runner_tmp --output spectra_pi0.csv \
    --z-points 16 --threads 4 --max-runners 4 \
    --pt-min 1 --pt-max 10 --pt-step 1 --clean
```
