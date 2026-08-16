    program polyroots_test_2

    !! Test polynomial
    !! 3*x**4 - 9*x**3 - 369*x**2 + 705*x + 4950
    !!
    !! Its exact roots are:
    !! -10, -3, 5, 11
    !!
    !! see: https://community.intel.com/t5/Intel-Fortran-Compiler/IFX-2026-1-returns-incorrect-RPOLY-roots-with-O1-and-O2-fp/td-p/1756601

    use polyroots_module, only: wp => polyroots_module_rk, rpoly

    implicit none

    real(wp), parameter :: coefficients(5) = [ &
        3.0_wp, -9.0_wp, -369.0_wp, 705.0_wp, 4950.0_wp]
    real(wp) :: real_parts(4), imaginary_parts(4)
    complex(wp) :: root, value
    real(wp) :: maximum_residual
    integer :: i, status, j

    call rpoly(coefficients, 4, real_parts, imaginary_parts, status)

    maximum_residual = 0.0_wp
    write(*, '(a,i0)') 'status: ', status
    do i = 1, 4
        root = cmplx(real_parts(i), imaginary_parts(i), wp)
        value = coefficients(1)
        do j = 2, 5
            value = value * root + coefficients(j) ! horner's rule
        end do
        maximum_residual = max(maximum_residual, abs(value))
        write(*, '(i0,2es25.16)') i, real_parts(i), imaginary_parts(i)
    end do
    write(*, '(a,es25.16)') 'maximum residual: ', maximum_residual

    if (status /= 0 .or. maximum_residual > 1.0e-8_wp) then
        error stop 'incorrect roots'
    end if

    end program polyroots_test_2