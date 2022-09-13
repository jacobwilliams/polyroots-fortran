![polyroots-fortran](media/logo.png)
============

**polyroots-fortran**: Polynomial Roots with Modern Fortran

[![GitHub release](https://img.shields.io/github/release/jacobwilliams/polyroots-fortran.svg)](https://github.com/jacobwilliams/polyroots-fortran/releases/latest)
[![CI Status](https://github.com/jacobwilliams/polyroots-fortran/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/polyroots-fortran/actions)
[![codecov](https://codecov.io/gh/jacobwilliams/polyroots-fortran/branch/master/graph/badge.svg)](https://codecov.io/gh/jacobwilliams/polyroots-fortran)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/polyroots-fortran)](https://github.com/jacobwilliams/polyroots-fortran/commits/master)

## Description

A modern Fortran library for finding the roots of polynomials.

**This is a work in progress**

## Methods

Method name | Polynomial type | Coefficients | Roots | Description
--- | --- | --- | --- | ---
[`rpoly`](https://jacobwilliams.github.io/polyroots-fortran/proc/rpoly.html) | General | real | complex | [Jenkins & Traub](https://dl.acm.org/doi/10.1145/355637.355643)
[`dcbcrt`](https://jacobwilliams.github.io/polyroots-fortran/proc/dcbcrt.html) | Cubic | real | complex | [NSWC Library](https://github.com/jacobwilliams/nswc)
[`dqdcrt`](https://jacobwilliams.github.io/polyroots-fortran/proc/dqdcrt.html) | Quadratic | real | complex | [NSWC Library](https://github.com/jacobwilliams/nswc)
[`qr_algeq_solver`](https://jacobwilliams.github.io/polyroots-fortran/proc/qr_algeq_solver.html) | Monic | real | complex | [Edelman & Murakami (1995)](https://www.ams.org/journals/mcom/1995-64-210/S0025-5718-1995-1262279-2/S0025-5718-1995-1262279-2.pdf)


## Compiling

A `fmp.toml` file is provided for compiling polyroots-fortran with the [Fortran Package Manager](https://github.com/fortran-lang/fpm). For example, to build:

```
fpm build --profile release
```

By default, the library is built with double precision (`real64`) real values. Explicitly specifying the real kind can be done using the following processor flags:

Preprocessor flag | Kind  | Number of bytes
----------------- | ----- | ---------------
`REAL32`  | `real(kind=real32)`  | 4
`REAL64`  | `real(kind=real64)`  | 8
`REAL128` | `real(kind=real128)` | 16

For example, to build a single precision version of the library, use:

```
fpm build --profile release --flag "-DREAL32"
```

To run the unit tests:

```
fpm test
```

To use `polyroots-fortran` within your fpm project, add the following to your `fpm.toml` file:
```toml
[dependencies]
polyroots-fortran = { git="https://github.com/jacobwilliams/polyroots-fortran.git" }
```

or, to use a specific version:
```toml
[dependencies]
polyroots-fortran = { git="https://github.com/jacobwilliams/polyroots-fortran.git", tag = "1.0.0"  }
```

To generate the documentation using [ford](https://github.com/Fortran-FOSS-Programmers/ford), run: ```ford polyroots-fortran.md```

## Documentation

The latest API documentation for the `master` branch can be found [here](https://jacobwilliams.github.io/polyroots-fortran/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

## License

The polyroots-fortran source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/polyroots-fortran/blob/master/LICENSE.md) (BSD-style).

## See also

* [Roots-Fortran](https://github.com/jacobwilliams/roots-fortran)

## Similar libraries in other programming languages

* R: [polyroot](https://stat.ethz.ch/R-manual/R-devel/library/base/html/polyroot.html)
* MATLAB: [roots](https://www.mathworks.com/help/matlab/ref/roots.html)
* C: [GSL - Polynomials](https://www.gnu.org/software/gsl/doc/html/poly.html), [MPSolve](https://numpi.dm.unipi.it/software/mpsolve)
* Julia: [PolynomialRoots.jl](https://github.com/giordano/PolynomialRoots.jl), [FastPolynomialRoots.jl](https://github.com/andreasnoack/FastPolynomialRoots.jl), [Polynomials.jl](https://github.com/JuliaMath/Polynomials.jl)
* Python: [numpy.polynomial.polynomial](https://docs.scipy.org/doc//numpy-1.10.4/reference/routines.polynomials.polynomial.html)


## References

* [GAMS Class F1a](https://gams.nist.gov/cgi-bin/serve.cgi/Class/F1a).
* `RPOLY` and `CPOLY` from Jenkins and Traub (1975) which can be found on Netlib. See also [Alan Miller's updated version](https://jblevins.org/mirror/amiller/rpoly.f90).
* `qralg.f` (part of [`/opt/companion.tgz`](https://netlib.org/opt/companion.tgz) on Netlib) from [Edelman & Murakami (1995)](https://www.ams.org/journals/mcom/1995-64-210/S0025-5718-1995-1262279-2/S0025-5718-1995-1262279-2.pdf),
* [Codes from Skowron & Gould (2012)](http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/) and
* [eiscor - eigensolvers based on unitary core transformations](https://github.com/eiscor/eiscor) containing the AMVW method from the work of [Aurentz et al. (2015), Fast and Backward Stable Computation of Roots of Polynomials](https://doi.org/10.1137/140983434) (an earlier version can be picked up from [the website of Ran Vandebril](https://people2.cs.kuleuven.be/~raf.vandebril/homepage/software/companion_qr.php?menu=5), one of the co-authors of that paper).
* [a solver using the companion matrix and LAPACK](https://fortran-lang.discourse.group/t/cardanos-solution-of-the-cubic-equation/111/5?u=ivanpribec)
* [Root-finding algorithms: Roots of Polynomials | Wikipedia](https://en.wikipedia.org/wiki/Root-finding_algorithms#Roots_of_polynomials)
* [Polynomial Roots | Wolfram MathWorld](https://mathworld.wolfram.com/PolynomialRoots.html)
* [What is a Companion Matrix | Nick Higham](https://nhigham.com/2021/03/23/what-is-a-companion-matrix/)
* [19 Dubious Ways to Compute the Zeros of a Polynomial | Cleve's Corner](https://blogs.mathworks.com/cleve/2016/06/27/19-dubious-ways-to-compute-the-zeros-of-a-polynomial/)
* [New Progress in Polynomial Root-finding | Victor Y. Pan](https://arxiv.org/pdf/1805.12042.pdf)

## See also

 * [Code coverage statistics](https://app.codecov.io/gh/jacobwilliams/polyroots-fortran) [codecov.io]
