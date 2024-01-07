![polyroots-fortran](media/logo.png)
============

**polyroots-fortran**: Polynomial Roots with Modern Fortran

[![GitHub release](https://img.shields.io/github/release/jacobwilliams/polyroots-fortran.svg)](https://github.com/jacobwilliams/polyroots-fortran/releases/latest)
[![CI Status](https://github.com/jacobwilliams/polyroots-fortran/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/polyroots-fortran/actions)
[![codecov](https://codecov.io/gh/jacobwilliams/polyroots-fortran/branch/master/graph/badge.svg)](https://codecov.io/gh/jacobwilliams/polyroots-fortran)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/polyroots-fortran)](https://github.com/jacobwilliams/polyroots-fortran/commits/master)

## Description

A modern Fortran library for finding the roots of polynomials.

## Methods

Many of the methods are from legacy libraries. They have been extensively modified and refactored into Modern Fortran.

Method name | Polynomial type | Coefficients | Roots | Reference
--- | --- | --- | --- | ---
[`cpoly`](https://jacobwilliams.github.io/polyroots-fortran/proc/cpoly.html) | General | complex | complex | [Jenkins & Traub (1972)](https://dl.acm.org/doi/abs/10.1145/361254.361262)
[`rpoly`](https://jacobwilliams.github.io/polyroots-fortran/proc/rpoly.html) | General | real | complex | [Jenkins & Traub (1975)](https://dl.acm.org/doi/10.1145/355637.355643)
[`rpzero`](https://jacobwilliams.github.io/polyroots-fortran/proc/rpzero.html) | General | real | complex | [SLATEC](https://netlib.org/slatec/src/rpzero.f)
[`cpzero`](https://jacobwilliams.github.io/polyroots-fortran/proc/cpzero.html) | General | complex | complex | [SLATEC](https://netlib.org/slatec/src/cpzero.f)
[`rpqr79`](https://jacobwilliams.github.io/polyroots-fortran/proc/rpqr79.html) | General | real | complex | [SLATEC](https://netlib.org/slatec/src/rpqr79.f)
[`cpqr79`](https://jacobwilliams.github.io/polyroots-fortran/proc/cpqr79.html) | General | complex | complex | [SLATEC](https://netlib.org/slatec/src/cpqr79.f)
[`dcbcrt`](https://jacobwilliams.github.io/polyroots-fortran/proc/dcbcrt.html) | Cubic | real | complex | [NSWC Library](https://github.com/jacobwilliams/nswc)
[`dqdcrt`](https://jacobwilliams.github.io/polyroots-fortran/proc/dqdcrt.html) | Quadratic | real | complex | [NSWC Library](https://github.com/jacobwilliams/nswc)
[`dpolz`](https://jacobwilliams.github.io/polyroots-fortran/proc/dpolz.html) | General | real | complex | [MATH77 Library](https://netlib.org/math/)
[`cpolz`](https://jacobwilliams.github.io/polyroots-fortran/proc/cpolz.html) | General | complex | complex | [MATH77 Library](https://netlib.org/math/)
[`polyroots`](https://jacobwilliams.github.io/polyroots-fortran/proc/polyroots.html) | General | real | complex | [LAPACK](https://netlib.org/lapack/explore-html/index.html)
[`cpolyroots`](https://jacobwilliams.github.io/polyroots-fortran/proc/cpolyroots.html) | General | complex | complex | [LAPACK](https://netlib.org/lapack/explore-html/index.html)
[`rroots_chebyshev_cubic`](https://jacobwilliams.github.io/polyroots-fortran/proc/rroots_chebyshev_cubic.html) | Cubic | real | complex | [Lebedev (1991)](https://doi.org/10.1515/rnam.1991.6.4.315)
[`qr_algeq_solver`](https://jacobwilliams.github.io/polyroots-fortran/proc/qr_algeq_solver.html) | General | real | complex | [Edelman & Murakami (1995)](https://www.ams.org/journals/mcom/1995-64-210/S0025-5718-1995-1262279-2/S0025-5718-1995-1262279-2.pdf)
[`polzeros`](https://jacobwilliams.github.io/polyroots-fortran/proc/polzeros.html) | General | complex | complex | [Bini (1996)](https://link.springer.com/article/10.1007/BF02207694)
[`cmplx_roots_gen`](https://jacobwilliams.github.io/polyroots-fortran/proc/cmplx_roots_gen.html) | General | complex | complex | [Skowron & Gould (2012)](http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/)
[`fpml`](https://jacobwilliams.github.io/polyroots-fortran/proc/fpml.html) | General | complex | complex | [Cameron (2019)](https://link.springer.com/article/10.1007/s11075-018-0641-9)

## Example

An example of using `polyroots` to compute the zeros for a 5th order real polynomial $$P(x) = x^5 + 2x^4 + 3x^3 + 4x^2 + 5x + 6$$

```fortran
program example

use iso_fortran_env
use polyroots_module, wp => polyroots_module_rk

implicit none

integer,parameter :: degree = 5 !! polynomial degree
real(wp),dimension(degree+1) :: p = [1,2,3,4,5,6] !! coefficients

integer :: i !! counter
integer :: istatus !! status code
real(wp),dimension(degree) :: zr !! real components of roots
real(wp),dimension(degree) :: zi !! imaginary components of roots

call polyroots(degree, p, zr, zi, istatus)

write(*,'(A,1x,I3)') 'istatus: ', istatus
write(*, '(*(a22,1x))') 'real part', 'imaginary part'
do i = 1, degree
    write(*,'(*(e22.15,1x))') zr(i), zi(i)
end do

end program example
```

The result is:

```
istatus:    0
             real part         imaginary part
 0.551685463458982E+00  0.125334886027721E+01
 0.551685463458982E+00 -0.125334886027721E+01
-0.149179798813990E+01  0.000000000000000E+00
-0.805786469389031E+00  0.122290471337441E+01
-0.805786469389031E+00 -0.122290471337441E+01
```

## Compiling

A `fpm.toml` file is provided for compiling polyroots-fortran with the [Fortran Package Manager](https://github.com/fortran-lang/fpm). For example, to build:

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

All routines, except for `polyroots` are available for any of the three real kinds. `polyroots` is not available for `real128` kinds since there is no corresponding LAPACK eigenvalue solver.

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
polyroots-fortran = { git="https://github.com/jacobwilliams/polyroots-fortran.git", tag = "1.2.0"  }
```

To generate the documentation using [ford](https://github.com/Fortran-FOSS-Programmers/ford), run: ```ford ford.md```

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

## Other references and codes

* [GAMS Class F1a](https://gams.nist.gov/cgi-bin/serve.cgi/Class/F1a).
* [eiscor - eigensolvers based on unitary core transformations](https://github.com/eiscor/eiscor) containing the AMVW method from the work of [Aurentz et al. (2015), Fast and Backward Stable Computation of Roots of Polynomials](https://doi.org/10.1137/140983434) (an earlier version can be picked up from [the website of Ran Vandebril](https://people2.cs.kuleuven.be/~raf.vandebril/homepage/software/companion_qr.php?menu=5), one of the co-authors of that paper).
* [PA03](https://www.hsl.rl.ac.uk/archive/specs/pa03.pdf) HSL Archive code for computing all the roots of a cubic polynomial
* [PA05](https://www.hsl.rl.ac.uk/archive/specs/pa05.pdf) HSL Archive code for computing all the roots of a quartic polynomial
* [PA16](https://www.hsl.rl.ac.uk/catalogue/pa16.html), [PA17](https://www.hsl.rl.ac.uk/catalogue/pa17.html) HSL codes for computing zeros of polynomials using method of Madsen and Reid
* Various codes from [Alan Miller](https://jblevins.org/mirror/amiller/)
* [A solver using the companion matrix and LAPACK](https://fortran-lang.discourse.group/t/cardanos-solution-of-the-cubic-equation/111/5?u=ivanpribec)
* [Root-finding algorithms: Roots of Polynomials | Wikipedia](https://en.wikipedia.org/wiki/Root-finding_algorithms#Roots_of_polynomials)
* [Polynomial Roots | Wolfram MathWorld](https://mathworld.wolfram.com/PolynomialRoots.html)
* [What is a Companion Matrix | Nick Higham](https://nhigham.com/2021/03/23/what-is-a-companion-matrix/)
* [19 Dubious Ways to Compute the Zeros of a Polynomial | Cleve's Corner](https://blogs.mathworks.com/cleve/2016/06/27/19-dubious-ways-to-compute-the-zeros-of-a-polynomial/)
* [New Progress in Polynomial Root-finding | Victor Y. Pan](https://arxiv.org/pdf/1805.12042.pdf)

## See also

 * [Code coverage statistics](https://app.codecov.io/gh/jacobwilliams/polyroots-fortran) [codecov.io]
