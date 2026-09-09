#pragma once

// Builds/tears down every per-kinematic-point interpolation table (the NLO
// coefficient functions, their xi-convolution, and the fixed-Y dipole
// slice) for the given (xp, xg, k). Must be called once before, and
// clear_interp() once after, each (z,k) point evaluated in main().
void init_interp(double xp, double xg, double k);
void clear_interp();
