!*****************************************************************************************
!>
!  Tests for [[rpoly]].

    program rpoly_test

    use iso_fortran_env
    use polyroots_module, wp => polyroots_module_rk

    implicit none

    integer,parameter :: max_degree = 10
    integer,parameter :: n_cases = 20 !! number of cases to run

    real(wp),dimension(max_degree+1) :: p, zr, zi, s
    complex(wp),dimension(max_degree+1) :: r, cp
    complex(wp),dimension(6*(max_degree+1)) :: t !! work array for [[rpzero]]
    real(wp) :: detil
    integer :: degree, i, istatus
    logical :: fail
    integer :: icase
    integer, allocatable :: seed(:)
    integer :: n

    !--------------------------------------

    ! set random seed for consistent results:
    call random_seed(size=n)
    allocate(seed(n))
    seed = 42
    call random_seed(put=seed)

    do icase = 1, n_cases

        write(*,'(/A,I2,A)') '--------CASE ', icase, ' ---------'

        if (icase==1) then
            degree = 10
            p(1)  = 1._wp
            p(2)  = -55._wp
            p(3)  = 1320._wp
            p(4)  = -18150._wp
            p(5)  = 157773._wp
            p(6)  = -902055._wp
            p(7)  = 3416930._wp
            p(8)  = -8409500._wp
            p(9)  = 12753576._wp
            p(10) = -10628640._wp
            p(11) = 3628800._wp
        else
            degree = get_random_integer_number(3,max_degree)
            p(1)  = get_random_number(-1000.0_wp,1000.0_wp)
            p(2)  = get_random_number(-1000.0_wp,1000.0_wp)
            p(3)  = get_random_number(-1000.0_wp,1000.0_wp)
            p(4)  = get_random_number(-1000.0_wp,1000.0_wp)
            p(5)  = get_random_number(-1000.0_wp,1000.0_wp)
            p(6)  = get_random_number(-1000.0_wp,1000.0_wp)
            p(7)  = get_random_number(-1000.0_wp,1000.0_wp)
            p(8)  = get_random_number(-1000.0_wp,1000.0_wp)
            p(9)  = get_random_number(-1000.0_wp,1000.0_wp)
            p(10) = get_random_number(-1000.0_wp,1000.0_wp)
            p(11) = get_random_number(-1000.0_wp,1000.0_wp)
        end if

        write(*,'(A,1X/,*(g23.15/))') ' Coefficients: ', p(1:degree+1)

        write(*, '(A,1x,i3)') 'rpoly'
        write(*, '(a)') '  real part               imaginary part         root'
        call rpoly(p, degree, zr, zi, fail)
        if (fail) error stop ' ** failure by rpoly **'
        call check_results(zr, zi, degree)

        write(*, '(/A,1x,i3)') 'rpzero'
        write(*, '(a)') '  real part               imaginary part         root'
        istatus = 0 ! no estimates input
        call rpzero(degree,p,r,t,istatus,s)
        if (istatus/=0) error stop ' ** failure by rpzero **'
        call check_results(real(r,wp), aimag(r), degree)

        write(*, '(/A,1x,i3)') 'rpqr79'
        write(*, '(a)') '  real part               imaginary part         root'
        call rpqr79(degree,p,r,istatus)
        if (istatus/=0) error stop ' ** failure by rpqr79 **'
        call check_results(real(r,wp), aimag(r), degree)

        ! for now, just test this one with the real coefficients only:
        write(*, '(/A,1x,i3)') 'cpqr79'
        write(*, '(a)') '  real part               imaginary part         root'
        do i = 1, degree+1
            cp(i) = cmplx(p(i), 0.0_wp, wp) ! put in a complex number
        end do
        call cpqr79(degree,cp,r,istatus)
        if (istatus/=0) error stop ' ** failure by cpqr79 **'
        call check_results(real(r,wp), aimag(r), degree)

        write(*, '(/A,1x,i3)') 'qr_algeq_solver'
        write(*, '(a)') '  real part               imaginary part         root'
        call qr_algeq_solver(degree,p,zr,zi,detil,istatus)
        call check_results(zr, zi, degree)

        if (wp /= REAL128) then
            write(*, '(/A,1x,i3)') 'polyroots'
            write(*, '(a)') '  real part               imaginary part         root'
            call polyroots(degree, p, zr, zi, istatus)
            if (istatus/=0) error stop 'error: polyroots did not converge'
            call check_results(zr, zi, degree)
        end if

    end do

    !--------------------------------------
    ! this test provided by larry wigton

    write(*, '(/A)') 'rpoly example 2. a case where 1 is an obvious root.'
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

    contains

    !*****************************************************************************************
        subroutine check_results(re, im, degree)

            !! check the results.
            !! if any are not within the tolerance,
            !! then also try to polish them using the newton method.

            real(wp),dimension(:),intent(in) :: re, im
            integer,intent(in) :: degree

            complex(wp) :: z, root
            integer :: i,j !! counter
            real(wp) :: zr, zi ! copy of inputs for polishing
            integer :: istat

            real(wp),parameter :: tol = 1.0e-2_wp  !! acceptable root tolerance for tests
            real(wp),parameter :: ftol = 1.0e-8_wp !! desired root tolerance
            real(wp),parameter :: ztol = 10*epsilon(1.0_wp) !! newton tol for x

            do j = 1, degree
                z = cmplx(re(j), im(j), wp)
                root = p(1)
                do i = 2, degree+1
                    root = root * z + p(i) ! horner's rule
                end do
                write(*, '(3(2g23.15,1x))') re(j), im(j), abs(root)
                if (abs(root) > ftol) then
                    ! attempt to polish the root:
                    zr = re(j)
                    zi = im(j)
                    call newton_root_polish(degree, p, zr, zi, &
                                            ftol=ftol, ztol=ztol, maxiter=10, &
                                            istat=istat)
                    z = cmplx(zr, zi, wp) ! recompute root with possibly updated values
                    root = p(1)
                    do i = 2, degree+1
                        root = root * z + p(i) ! horner's rule
                    end do
                    write(*, '(3(2g23.15,1x),1X,A)') zr, zi, abs(root), 'POLISHED'
                    if (abs(root) > tol) error stop 'Error: insufficient accuracy'
                end if
            end do

        end subroutine check_results
    !*****************************************************************************************

    !*****************************************************************************************
    !> author: Jacob Williams
    !
    !  Returns a uniform random number `x`, such that: `a <= x < b`.
    !
    !  This routine is from the Fortran Astrodynamics Toolkit.

        function get_random_number(a,b) result(x)

            implicit none

            real(wp)            :: x
            real(wp),intent(in) :: a
            real(wp),intent(in) :: b

            call random_number(x)

            x = a + (b-a)*x

        end function get_random_number
    !*****************************************************************************************

    !*****************************************************************************************
    !>
    !  Returns a uniform integer random number `x`, such that: `a <= x < b`.

        function get_random_integer_number(a,b) result(i)

            implicit none

            integer            :: i
            integer,intent(in) :: a
            integer,intent(in) :: b

            real(wp) :: x

            call random_number(x)

            i = ceiling(a + (b-a)*x)

        end function get_random_integer_number
    !*****************************************************************************************

!*****************************************************************************************
    end program rpoly_test
!*****************************************************************************************