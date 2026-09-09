#pragma once

// Tiny shared helper used throughout the codebase. Split out on its own so
// every translation unit that needs it can pull in a single trivial header
// instead of a heavier one.
static inline double Sq(double x){ return x*x; }
