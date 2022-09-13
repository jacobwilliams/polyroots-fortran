!*****************************************************************************************
!>
!  Tests for [[rpoly]].

    program rpoly_test

    use iso_fortran_env
    use polyroots_module, only: wp => polyroots_module_rk, rpoly

    implicit none

    real (wp)  :: p(50), zr(50), zi(50)
    integer    :: degree, i
    logical    :: fail

    write(*, '(A)') ' example 1. polynomial with zeros 1,2,...,10.'

    degree = 10
    p(1) = 1._wp
    p(2) = -55._wp
    p(3) = 1320._wp
    p(4) = -18150._wp
    p(5) = 157773._wp
    p(6) = -902055._wp
    p(7) = 3416930._wp
    p(8) = -8409500._wp
    p(9) = 12753576._wp
    p(10) = -10628640._wp
    p(11) = 3628800._wp

    call rpoly(p, degree, zr, zi, fail)
    if (fail) then
        error stop ' ** failure by rpoly **'
    else
        write(*, '(a/ (2g23.15))') ' real part           imaginary part',  &
                                    (zr(i), zi(i), i=1,degree)
    end if

    ! this test provided by larry wigton

    write(*, *)
    write(*, *) "now try case where 1 is an obvious root"
    degree = 5
    p(1) = 8.d0
    p(2) = -8.d0
    p(3) = 16.d0
    p(4) = -16.d0
    p(5) = 8.d0
    p(6) = -8.d0

    call rpoly(p, degree, zr, zi, fail)
    if (fail) then
        error stop ' ** failure by rpoly **'
    else
        write(*, *) ' real part           imaginary part'
        write(*, '(2g23.15)') (zr(i), zi(i), i=1,degree)
    end if

!*****************************************************************************************
    end program rpoly_test
!*****************************************************************************************