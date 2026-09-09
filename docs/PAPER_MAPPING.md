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

## The coefficient functions: name mismatches

The nine functions in [nlo_coefficients.hpp](../src/nlo_coefficients.hpp)
implement the transverse-integral terms of Eqs. (10) and (12). Most keep
the paper's own symbol, but **three don't** — this is the single most
important thing to know when cross-referencing code and paper:

| Paper symbol | Equation | Code symbol | Relation |
|---|---|---|---|
| I₁ | (10c) | `I1` (free function, [point_tables.cpp](../src/point_tables.cpp)) | same |
| I₂ | (10d) | `NLOCoefficients::I2` | same |
| J | (10a) | `NLOCoefficients::J` | same |
| Jᵥ | (10b) | `NLOCoefficients::Jv` | same |
| H₁ | (12a) | `H1` (free function, [point_tables.cpp](../src/point_tables.cpp)) | same |
| H₂ | (12b) | `NLOCoefficients::H2` | same |
| H₃ | (12c) | `NLOCoefficients::H3` | same |
| **H₄** | (12d) | **`NLOCoefficients::H5`** | same function, different name |
| **K₁** | (12e) | **`NLOCoefficients::J1`** | `J1(r,ξ) == (1/4)·K₁(r,ξ)` |
| **K₂** | (12f) | **`NLOCoefficients::K3`** | `K3(r,ξ) == (1/4)·K₂(r,ξ)` |

I.e. the code's `J1` has nothing to do with the paper's `J`, and the
code's `K3` is unrelated to any "K3" in the paper (there is no K3 in the
paper — only K1 and K2). Best guess at how this happened: an earlier,
larger set of candidate terms was probably tried during development (the
`flag` values 0-7 in `nlo_coefficients.cpp`'s `integrand_x` hint at more
variants than ended up in the final paper), and the surviving ones kept
their original working names rather than being renamed to match the
paper on publication.

`I1`/`H1` need no transverse integral (they're plain algebra in S(r,Y)),
so they aren't part of the `NLOCoefficients` class — they're free
functions in `point_tables.cpp`, right next to the channel decomposition
that uses them.

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
