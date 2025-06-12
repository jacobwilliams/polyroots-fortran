!*****************************************************************************************
!>
!  Example in the readme.

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

! if we have a quad solver:
#ifdef USE_STDLIB_LAPACK

call polyroots(degree, p, zr, zi, istatus)

write(*,'(/A,1x,I3)') 'istatus: ', istatus
write(*, '(*(a22,1x))') 'real part', 'imaginary part', 'root'
do i = 1, degree
    write(*,'(*(e22.15,1x))') zr(i), zi(i), abs(evaluate(p, cmplx(zr(i), zi(i), wp)))
end do

#endif

contains

    function evaluate(p, x) result(f)

        !! evaluate the polynomial:
        !! `f = p(1)*x**n + p(2)*x**(n-1) + ... + p(n)*x + p(n+1)`
        !! and its derivative using Horner's Rule.
        !!
        !! See: "Roundoff in polynomial evaluation", W. Kahan, 1986
        !! https://rosettacode.org/wiki/Horner%27s_rule_for_polynomial_evaluation#Fortran

        implicit none

        real(wp),dimension(:),intent(in) :: p !! coefficients
        complex(wp), intent(in) :: x
        complex(wp) :: f    !! function value at `x`

        integer :: i !! counter

        f = p(1)
        do i = 2, size(p)-1 + 1
            f = f*x + p(i)
        end do

    end function evaluate

end program example
!*****************************************************************************************