!*****************************************************************************************
!>
!  Tests for [[rpoly]].

    program rpoly_test

    use iso_fortran_env
    use polyroots_module, wp => polyroots_module_rk

    implicit none

    real (wp)  :: p(50), zr(50), zi(50), detil, s(50)
    integer    :: degree, i, istatus
    logical    :: fail
    complex(wp) :: r(50)
    complex(wp) :: t(66) !! work array for rpzero

    !--------------------------------------

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

    write(*, '(/A)') 'rpoly example 1. polynomial with zeros 1,2,...,10.'
    call rpoly(p, degree, zr, zi, fail)
    if (fail) error stop ' ** failure by rpoly **'
    write(*, '(a/ (2g23.15))') ' real part           imaginary part',  &
                                    (zr(i), zi(i), i=1,degree)

    write(*, '(/A)') 'rpzero example 1. polynomial with zeros 1,2,...,10.'
    istatus = 0 ! no estimates input
    call rpzero(degree,p,r,t,istatus,s)
    if (istatus/=0) error stop ' ** failure by rpzero **'
    write(*, '(a/ (2g23.15))') ' real part           imaginary part',  &
                                    (real(r(i),wp), aimag(r(i)), i=1,degree)

    write(*, '(/A)') 'rpqr79 example 1. polynomial with zeros 1,2,...,10.'
    call rpqr79(degree,p,r,istatus)
    if (istatus/=0) error stop ' ** failure by rpqr79 **'
    write(*, '(a/ (2g23.15))') ' real part           imaginary part',  &
                                    (real(r(i),wp), aimag(r(i)), i=1,degree)

    write(*, '(/A)') 'qr_algeq_solver example 1. polynomial with zeros 1,2,...,10.'
    call qr_algeq_solver(degree,p,zr,zi,detil,istatus)
    write(*, '(a/ (2g23.15))') ' real part           imaginary part',  &
            (zr(i), zi(i), i=1,degree)

    write(*, '(/A)') 'polyroots example 1. polynomial with zeros 1,2,...,10.'
    call polyroots(degree, p, zr, zi, istatus)
    if (istatus/=0) error stop 'error: polyroots did not converge'
    write(*, '(a/ (2g23.15))') ' real part           imaginary part',  &
            (zr(i), zi(i), i=1,degree)

    !--------------------------------------
    ! this test provided by larry wigton

    write(*, '(/A)') 'rpoly example 2. polynomial with zeros 1,2,...,10.'
    write(*, '(A)')  '[a case where 1 is an obvious root]'
    degree = 5
    p(1) = 8.0_wp
    p(2) = -8.0_wp
    p(3) = 16.0_wp
    p(4) = -16.0_wp
    p(5) = 8.0_wp
    p(6) = -8.0_wp

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