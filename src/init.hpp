#pragma once

// One-time setup: loads the LHAPDF set, reads and interpolates the BK
// solution, and installs the GSL error handler. Must be called once before
// any per-kinematic-point computation (init_interp / sigma_LO_k /
// sigma_NLO_k).
void init();
