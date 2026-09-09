# Paper ↔ code mapping

This codebase implements the calculation of:

> Heikki Mäntysaari and Yossathorn Tawabutr, *"Complete Next-to-Leading Order
> Calculation of Single Inclusive π⁰ Production in Forward Proton-Nucleus
> Collisions"*, [arXiv:2310.06640](https://arxiv.org/abs/2310.06640).

This document maps the paper's equations and notation onto the code's
classes and functions. It was produced by extracting the paper's text
(`pdftotext -layout`) and checking each formula below against the
corresponding code term-by-term — see "Confidence" per section. If you
spot a mismatch, trust the paper over this document and fix this file.

## Section map

| Paper section | Content | Code |
|---|---|---|
| Sec. II, Eq. (1) | Factorization: dσ = f ⊗ H ⊗ S ⊗ D | Not fully implemented: no fragmentation function D (hadron level, sec. V) — see "Parton vs. hadron level" below |
| Sec. II, Eqs. (4)-(6) | Kinematics: X(ξ), Xg, xp | [main.cpp](../src/main.cpp) (xp, xg), [point_tables.cpp](../src/point_tables.cpp)'s `integrand_xi` (X(ξ), y) |
| Sec. II, Eqs. (7a)/(7b) | LO cross section, quark/gluon channel | [sigma_LO.cpp](../src/sigma_LO.cpp)'s `sigma_LO_r` |
| Sec. II, Eqs. (9)-(10) | NLO qq channel + I1/I2/J/Jv | [point_tables.cpp](../src/point_tables.cpp)'s `integrand_xi` (`with_CF`/`with_Nc`), [nlo_coefficients.cpp](../src/nlo_coefficients.cpp) |
| Sec. II, Eqs. (11)-(12) | NLO gg/qg/gq channels + H1-H4/K1/K2 | same, `with_gg`/`with_gl`/`with_gq` |
| Sec. III, Eqs. (13)-(14) | Dipole amplitude initial condition (pp / pA) | [dipole_amplitude.cpp](../src/dipole_amplitude.cpp)'s `Sr_0` |
| Sec. III | BK-evolved S(r,Y), read from a solved-BK grid | [dipole_amplitude.cpp](../src/dipole_amplitude.cpp)'s `DipoleAmplitude` |
| Sec. IV | Parton-level results (no FF convolution) | This is what `main.cpp`'s loop actually computes |
| Sec. V | Hadron-level results (with FF convolution) | **Not implemented** — see below |
| Appendix A | Running-coupling scheme comparison | [running_coupling.cpp](../src/running_coupling.cpp), selected via `RunParameters::alpha_s_running` |

## The coefficient functions

The nine functions in [nlo_coefficients.hpp](../src/nlo_coefficients.hpp)
implement the transverse-integral terms of Eqs. (10) and (12), and all
now use the paper's own symbols:

| Paper symbol | Equation | Code symbol |
|---|---|---|
| I₁ | (10c) | `I1` (free function, [point_tables.cpp](../src/point_tables.cpp)) |
| I₂ | (10d) | `NLOCoefficients::I2` |
| J | (10a) | `NLOCoefficients::J` |
| Jᵥ | (10b) | `NLOCoefficients::Jv` |
| H₁ | (12a) | `H1` (free function, [point_tables.cpp](../src/point_tables.cpp)) |
| H₂ | (12b) | `NLOCoefficients::H2` |
| H₃ | (12c) | `NLOCoefficients::H3` |
| H₄ | (12d) | `NLOCoefficients::H4` |
| K₁ | (12e) | `NLOCoefficients::K1` |
| K₂ | (12f) | `NLOCoefficients::K2` |

`I1`/`H1` need no transverse integral (they're plain algebra in S(r,Y)),
so they aren't part of the `NLOCoefficients` class — they're free
functions in `point_tables.cpp`, right next to the channel decomposition
that uses them.

### Former `J1`/`K3`/`H5` naming (historical)

Before this table was written, these three were named `J1`, `K3` and
`H5` — an unrelated, earlier naming with no connection to the paper's own
K1/K2/H4 numbering (`J1` had nothing to do with the paper's `J`, and
there was no "K3" in the paper at all — only K1 and K2). Best guess at
how that happened: an earlier, larger set of candidate terms was
probably tried during development (the `flag` values 0-7 in
`nlo_coefficients.cpp`'s `integrand_x` hint at more variants than ended
up in the final paper), and the surviving ones kept their original
working names rather than being renamed to match the paper on
publication.

Renaming `H5`→`H4` was a pure rename (same value). `J1`→`K1` and
`K3`→`K2` needed rescaling as well as renaming, since the old functions
computed exactly `1/4` of the paper's K1/K2 (their callers folded the
missing factor of 4 into whatever multiplied them, e.g.
`resgl2 = -tables.J1(r,y)` implicitly encoded the paper's `-(1/4)K1(r,ξ)`
term with the `1/4` hidden inside `J1` itself): `K1`/`K2` now return the
paper's literal Eq. (12e)/(12f) values, and the `1/4` was moved out to
the call sites in `point_tables.cpp`'s `integrand_xi`
(`resgl2 = -0.25*tables.K1(r,y)`, `resgq2 = -0.25*tables.K2(r,y)`),
matching Eqs. (11b)/(11c) as literally written. Purely a rename +
where-the-constant-lives change — every intermediate and final numeric
value is unchanged; verified against `tests/regression/`.

## Channel → equation → `RunParameters` flag

| Channel | Equation | Flag (set in `params::make_run_parameters`) |
|---|---|---|
| qq | (9) | `with_CF` and `with_Nc` (both true together — see below) |
| gg | (11a) | `with_gg` |
| qg | (11b) | `with_gl` |
| gq | (11c) | `with_gq` |

Eq. (9) is written in the large-N_c limit and its terms have different
color-factor scaling: the I₁/I₂ terms scale with C_F, the J/Jᵥ terms with
N_c. The code computes these as two separate contributions
(`with_CF`/`with_Nc`), both enabled together whenever the incoming and
outgoing partons are both quarks (see `make_run_parameters` in
[params.cpp](../src/params.cpp)) — this is a decomposition of Eq. (9), not
two different channels.

## Not part of Eqs. (9)-(12): `Jv2` and `JJv_xi1`

`NLOCoefficients::Jv2` and `NLOCoefficients::JJv_xi1` don't correspond to
anything in this paper's equations. `JJv_xi1` is gated behind
`params::with_xi1`, a compile-time `false`, so it's never actually
evaluated in this codebase. Based on the surrounding comments ("subtracted
scheme" vs. this paper's "unsubtracted scheme", sec. II) and the paper's
citation of it, these are most likely leftover machinery from the
alternative subtraction scheme of:

> B. Ducloué, E. Iancu, T. Lappi, A. H. Mueller, G. Soyez,
> D. N. Triantafyllopoulos and Y. Zhu, *"Use of a running coupling in the
> NLO calculation of forward hadron production"*, Phys. Rev. D 97 (2018)
> 054020, [arXiv:1712.07480](https://arxiv.org/abs/1712.07480)
> (ref. [26] in arXiv:2310.06640).

**This has not been verified** against arXiv:1712.07480's own equations —
flagging it here so it isn't mistaken for settled fact.

## Known coupling-scheme gaps (enforced at runtime)

Checking each channel's `alpha_s_running` handling against every one of
the seven schemes (`fixed`/`mom`/`parent`/`daughter`/`smallest`/`mixed`/
`mixedbd`) turned up three coefficient functions whose per-scheme
branches don't cover every scheme, where the gap is *not* filled by the
calling code either — meaning the affected term is silently computed
with no running-coupling factor at all (as if α_s = 1) rather than
failing loudly:

| Function | Channel | Covered internally ([nlo_coefficients.cpp](../src/nlo_coefficients.cpp)) | Covered externally ([point_tables.cpp](../src/point_tables.cpp)) | Gap |
|---|---|---|---|---|
| `K1` | qg (`with_gl`) | none — the alpha_s code for this branch is entirely commented out | `parent`, `mixed`, `mixedbd` | **`daughter`, `smallest`** |
| `K2` | gq (`with_gq`) | `daughter`, `mixed`, `parent` | `mixed`, `mixedbd` | **`smallest`** |
| `H2` | gg (`with_gg`) | `daughter`, `parent`, `mixed` | `mixedbd` | **`smallest`** |

The qq channel (`I1`/`I2`/`J`/`Jv`, `with_CF`/`with_Nc`) was checked too
and has no gap — every scheme is handled in all four functions.

`params::make_run_parameters` now calls a `validate_alpha_s_running`
check and exits with an error (rather than silently producing physically
wrong results) for:
- `alpha_s_running=smallest` combined with any of the qg/gq/gg channels
  (`incoming`/`outgoing` selecting `with_gl`/`with_gq`/`with_gg`), and
- `alpha_s_running=daughter` combined with the qg channel (`with_gl`).

None of this affects anything actually run so far: both `Script.sh` and
`Script_oberon.sh` hardcode `rc=mom`, and `mom` (like `fixed`) bypasses
every per-term branch entirely — it multiplies the whole assembled
xi-convolution by a single constant coupling in `sigma_NLO_r`, so it was
never exposed to this gap. `parent` was also checked and is complete for
all four channels.

## Possible missing σ0/2 normalization for `pp`

Eq. (8) states `∫d²b⊥ → σ0/2` for proton targets — a constant
multiplicative factor that should appear in every `pp` cross section
(Eqs. 7a/7b, 9, 11a-c). Tracing the code's normalization end-to-end (the
Hankel-transform convention in `sigma_LO_k`/`sigma_NLO_k` reproduces the
paper's `1/(4π²)` prefactor exactly), `params::sigma0` is referenced in
exactly one place in the entire codebase:

```
dipole_amplitude.cpp: Sr_0's pA branch (Eq. 14), via sigma0*Anucleus*rp.TA
```

It is **not** referenced anywhere in the `pp` path (`Sr_0`'s pp branch,
Eq. 13, correctly has no σ0 either — matching the paper). So the pA cross
section correctly carries its own σ0/2 (via Eq. 14's `A·T_A(b)`
construction), while the pp cross section appears to be missing the
Eq. (8) σ0/2 factor entirely, making every pp result too small by a fixed
constant (~47, for the currently-active σ0=94.46).

**This has deliberately not been "fixed" or guarded against** — unlike
the coupling-scheme gaps above, it isn't a (setting, channel) combination
that can be rejected; it's a question of whether σ0/2 is applied
elsewhere (a downstream plotting/analysis script) that only whoever
maintains those scripts can answer. If you confirm it's a genuine gap,
the fix is a single multiplication by `sigma0/2` somewhere in the pp path
(e.g. in `sigma_LO_r`/`sigma_NLO_r`, guarded on `rp.col != "pA"`).

## Parton vs. hadron level

The paper computes both parton-level (sec. IV) and hadron-level (sec. V,
via a DSS fragmentation-function convolution over the momentum fraction
z) cross sections. This codebase only computes the parton-level
quantities: `main.cpp`'s CLI takes `z` as a direct scan parameter
(`k = p/z`) rather than integrating over it, and there is no
fragmentation-function code at all (an earlier, fully-commented-out
attempt at this — calling into a DSS Fortran routine `fdss_` — was
removed as dead code during the RAII refactor; see the git history if you
need to resurrect it as a starting point).

## Known pre-existing issues (unrelated to this mapping, kept here for
## visibility)

- `col == "pA"` is currently broken: see
  [tests/regression/README.md](../tests/regression/README.md).
- A handful of latent warts flagged along the way (an uninitialized
  variable on unrecognized PDF flavor strings, a couple of unused
  variables, an unused parameter) are listed in the commit messages for
  the refactor steps that found them; none affect the results for any
  currently-used configuration.
