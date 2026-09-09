#pragma once

/*
DE-Quadrature
Numerical Automatic Integrator for Improper Integral
    method    : Double Exponential (DE) Transformation
    dimension : one
    table     : not use
functions
    intde  : integrator of f(x) over (a,b).
    intdei : integrator of f(x) over (a,infinity),
                 f(x) is non oscillatory function.
    intdeo : integrator of f(x) over (a,infinity),
                 f(x) is oscillatory function.

Third-party numerics (Ooura's DE-quadrature), unrelated to the physics in
the rest of this codebase -- kept as-is algorithmically. The only change
from the original is that f now also receives a `void *userdata` context
pointer, mirroring GSL's own gsl_function convention (already used
elsewhere in this codebase) instead of forcing every caller to smuggle
its inputs through file-scope globals. Pass nullptr if the integrand
needs no context.
*/

/*
intde
    [description]
        I = integral of f(x) over (a,b)
    [declaration]
        void intde(double (*f)(double, void*), double a, double b, double eps,
            double *i, double *err, void *userdata);
    [usage]
        intde(f, a, b, eps, &i, &err, userdata);
    [parameters]
        f         : integrand f(x, userdata) (double (*f)(double, void*))
        a         : lower limit of integration (double)
        b         : upper limit of integration (double)
        eps       : relative error requested (double)
        i         : approximation to the integral (double *)
        err       : estimate of the absolute error (double *)
        userdata  : context pointer passed through to f unchanged (void *)
    [remarks]
        function
            f(x) needs to be analytic over (a,b).
        relative error
            eps is relative error requested excluding
            cancellation of significant digits.
            i.e. eps means : (absolute error) /
                             (integral_a^b |f(x)| dx).
            eps does not mean : (absolute error) / I.
        error message
            err >= 0 : normal termination.
            err < 0  : abnormal termination (m >= mmax).
                       i.e. convergent error is detected :
                           1. f(x) or (d/dx)^n f(x) has
                              discontinuous points or sharp
                              peaks over (a,b).
                              you must divide the interval
                              (a,b) at this points.
                           2. relative error of f(x) is
                              greater than eps.
                           3. f(x) has oscillatory factor
                              and frequency of the oscillation
                              is very high.
*/
void intde(double (*f)(double, void*), double a, double b, double eps,
    double *i, double *err, void *userdata = nullptr);

/*
intdei
    [description]
        I = integral of f(x) over (a,infinity),
            f(x) has not oscillatory factor.
    [declaration]
        void intdei(double (*f)(double, void*), double a, double eps,
            double *i, double *err, void *userdata);
    [usage]
        intdei(f, a, eps, &i, &err, userdata);
    [parameters]
        f         : integrand f(x, userdata) (double (*f)(double, void*))
        a         : lower limit of integration (double)
        eps       : relative error requested (double)
        i         : approximation to the integral (double *)
        err       : estimate of the absolute error (double *)
        userdata  : context pointer passed through to f unchanged (void *)
    [remarks]
        function
            f(x) needs to be analytic over (a,infinity).
        relative error
            eps is relative error requested excluding
            cancellation of significant digits.
            i.e. eps means : (absolute error) /
                             (integral_a^infinity |f(x)| dx).
            eps does not mean : (absolute error) / I.
        error message
            err >= 0 : normal termination.
            err < 0  : abnormal termination (m >= mmax).
                       i.e. convergent error is detected :
                           1. f(x) or (d/dx)^n f(x) has
                              discontinuous points or sharp
                              peaks over (a,infinity).
                              you must divide the interval
                              (a,infinity) at this points.
                           2. relative error of f(x) is
                              greater than eps.
                           3. f(x) has oscillatory factor
                              and decay of f(x) is very slow
                              as x -> infinity.
*/
void intdei(double (*f)(double, void*), double a, double eps,
    double *i, double *err, void *userdata = nullptr);

/*
intdeo
    [description]
        I = integral of f(x) over (a,infinity),
            f(x) has oscillatory factor :
            f(x) = g(x) * sin(omega * x + theta) as x -> infinity.
    [declaration]
        void intdeo(double (*f)(double, void*), double a, double omega,
            double eps, double *i, double *err, void *userdata);
    [usage]
        intdeo(f, a, omega, eps, &i, &err, userdata);
    [parameters]
        f         : integrand f(x, userdata) (double (*f)(double, void*))
        a         : lower limit of integration (double)
        omega     : frequency of oscillation (double)
        eps       : relative error requested (double)
        i         : approximation to the integral (double *)
        err       : estimate of the absolute error (double *)
        userdata  : context pointer passed through to f unchanged (void *)
    [remarks]
        function
            f(x) needs to be analytic over (a,infinity).
        relative error
            eps is relative error requested excluding
            cancellation of significant digits.
            i.e. eps means : (absolute error) /
                             (integral_a^R |f(x)| dx).
            eps does not mean : (absolute error) / I.
        error message
            err >= 0 : normal termination.
            err < 0  : abnormal termination (m >= mmax).
                       i.e. convergent error is detected :
                           1. f(x) or (d/dx)^n f(x) has
                              discontinuous points or sharp
                              peaks over (a,infinity).
                           2. relative error of f(x) is
                              greater than eps.
*/
void intdeo(double (*f)(double, void*), double a, double omega, double eps,
    double *i, double *err, void *userdata = nullptr);
