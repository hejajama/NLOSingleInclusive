#pragma once

// GSL error handler installed in main(): ignores a few benign GSL warnings
// (max-iterations/tolerance/underflow/rounding) and aborts on anything else.
void gsl_error_handler(const char * reason,
                       const char * file,
                       int line,
                       int error);
