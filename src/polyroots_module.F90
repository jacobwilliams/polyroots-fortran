!*****************************************************************************************
!>
!  Polynomial Roots with Modern Fortran
!
!### Author
!  * Jacob Williams
!
!@note The default real kind (`wp`) can be
!      changed using optional preprocessor flags.
!      This library was built with real kind:
#ifdef REAL32
!      `real(kind=real32)` [4 bytes]
#elif REAL64
!      `real(kind=real64)` [8 bytes]
#elif REAL128
!      `real(kind=real128)` [16 bytes]
#else
!      `real(kind=real64)` [8 bytes]
#endif

    module polyroots_module

    use iso_fortran_env

    implicit none

    private

#ifdef REAL32
    integer,parameter,public :: polyroots_module_rk = real32   !! real kind used by this module [4 bytes]
#elif REAL64
    integer,parameter,public :: polyroots_module_rk = real64   !! real kind used by this module [8 bytes]
#elif REAL128
    integer,parameter,public :: polyroots_module_rk = real128  !! real kind used by this module [16 bytes]
#else
    integer,parameter,public :: polyroots_module_rk = real64   !! real kind used by this module [8 bytes]
#endif

    integer,parameter :: wp = polyroots_module_rk  !! local copy of `polyroots_module_rk` with a shorter name


    contains
!*****************************************************************************************


!*****************************************************************************************
    end module polyroots_module
!*****************************************************************************************