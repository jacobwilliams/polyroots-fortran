![polyroots-fortran](media/logo.png)
============

**polyroots-fortran**: Polynomial Roots with Modern Fortran

[![GitHub release](https://img.shields.io/github/release/jacobwilliams/polyroots-fortran.svg)](https://github.com/jacobwilliams/polyroots-fortran/releases/latest)
[![CI Status](https://github.com/jacobwilliams/polyroots-fortran/actions/workflows/CI.yml/badge.svg)](https://github.com/jacobwilliams/polyroots-fortran/actions)
[![codecov](https://codecov.io/gh/jacobwilliams/polyroots-fortran/branch/master/graph/badge.svg)](https://codecov.io/gh/jacobwilliams/polyroots-fortran)
[![last-commit](https://img.shields.io/github/last-commit/jacobwilliams/polyroots-fortran)](https://github.com/jacobwilliams/polyroots-fortran/commits/master)

## Description

A modern Fortran library for finding the roots of polynomials.

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

## Usage

## Documentation

The latest API documentation for the `master` branch can be found [here](https://jacobwilliams.github.io/polyroots-fortran/). This was generated from the source code using [FORD](https://github.com/Fortran-FOSS-Programmers/ford).

## License

The polyroots-fortran source code and related files and documentation are distributed under a permissive free software [license](https://github.com/jacobwilliams/polyroots-fortran/blob/master/LICENSE.md) (BSD-style).

## Similar libraries in other programming languages


## References

## See also

 * [Code coverage statistics](https://app.codecov.io/gh/jacobwilliams/polyroots-fortran) [codecov.io]
