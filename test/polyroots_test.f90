!*****************************************************************************************
!>
!  Tests for [[rpoly]].

    program polyroots_test

    use iso_fortran_env
    use polyroots_module, wp => polyroots_module_rk
    use mt19937_64
    !use eiscor_module, only: z_poly_roots

    implicit none

    integer,parameter :: n_cases = 14 + 110 !! number of cases to run

    real(wp),dimension(:),allocatable :: p, pi, zr, zi, s, q, radius, berr,cond
    integer,dimension(:),allocatable :: conv
    complex(wp),dimension(:),allocatable :: r, cp, cp_
    integer :: degree, i, istatus, icase, n
    !integer,dimension(:),allocatable :: seed
    real(wp) :: detil
    logical :: fail
    logical :: failure !! if any of the tests failed
    logical,dimension(:),allocatable :: err
    integer :: idegree !! counter for degrees to test
    integer :: n_degree !! number of tests run for each degree so far
    type(mt19937) :: rand !! for random number generation

    failure = .false.

    ! set random seed for consistent results:
    call rand%initialize(42)
    ! call random_seed(size=n)
    ! allocate(seed(n))
    ! seed = 42
    ! call random_seed(put=seed)
    idegree = 0
    n_degree = 1

    do icase = 1, n_cases

        write(*,'(/A,I3,A)') '--------CASE ', icase, ' ---------'

        select case (icase)
        case(1)
            call allocate_arrays(10)
            p  = [1._wp, &
                 -55._wp, &
                 1320._wp, &
                 -18150._wp, &
                 157773._wp, &
                 -902055._wp, &
                 3416930._wp, &
                 -8409500._wp, &
                 12753576._wp, &
                 -10628640._wp, &
                 3628800._wp ]
            pi = 0.0_wp
        case(2)
            call allocate_arrays(4)
            p = [1,-3,20,44,54]
            pi = 0.0_wp
        case(3)
            call allocate_arrays(6)
            p = [1,-2,2,1,6,-6,8]
            pi = 0.0_wp
        case(4)
            call allocate_arrays(5)
            p = [1,1,-8,-16,7,15]
            pi = 0.0_wp
        case(5)
            call allocate_arrays(5)
            p = [1,7,5,6,3,2]
            pi = 0.0_wp
        case(6)
            call allocate_arrays(5)
            p = [2,3,6,5,7,1]
            pi = 0.0_wp
        case(7)
            call allocate_arrays(6)
            p = [1,0,-14,0,49,0,-36]
            pi = 0.0_wp
        case(8)
            call allocate_arrays(8)
            p = [1,0,-30,0,273,0,-820,0,576]
            pi = 0.0_wp
        case(9)
            call allocate_arrays(4)
            p = [1,0,0,0,-16]
            pi = 0.0_wp
        case(10)
            call allocate_arrays(6)
            p = [1,-2,2,1,6,-6,8]
            pi = 0.0_wp
        case(11)
            ! a case where 1 is an obvious root
            call allocate_arrays(5)
            pi = 0.0_wp
            p = [8,-8,16,-16,8,-8]
        case(12)
            call allocate_arrays(3)
            p = [ -8.0e18_wp,3.0e12_wp,5.0e6_wp,1.0_wp]
            pi = 0.0_wp
        case(13)
            call allocate_arrays(3)
            p = [4.0_wp, 3.0_wp, 2.0_wp, 1.0_wp]
            pi = 0.0_wp
        case(14)
            call allocate_arrays(2)
            p = [3.0_wp, 2.0_wp, 1.0_wp]
            pi = 0.0_wp

        ! case(15)   ! case 90 when compiled with ifort

        !     ! produces a root that doesn't evaluate to zero
        !     ! (same result from numpy)
        !     call allocate_arrays(7)
        !     p = [ 6.60460235615585_wp,&
        !           935.171169456812_wp,&
        !           867.901578904887_wp,&
        !           352.381787706374_wp,&
        !           320.264380528809_wp,&
        !          -332.592394794503_wp,&
        !          -398.892469194654_wp,&
        !           22.9384139877562_wp ]
        !     pi = 0.0_wp

        case default
            ! test a set of random coefficients for each degree:
            if (idegree>10) then
                idegree = 0
                n_degree = n_degree + 1
            end if
            idegree = idegree + 1
            call allocate_arrays(n_degree)

            do i = 1, degree+1
                p(i)  = get_random_number(-1000.0_wp,1000.0_wp)
                pi(i) = get_random_number(-10000.0_wp,10000.0_wp)
            end do
        end select
        do i = 1, degree+1
            cp(i) = cmplx(p(i), pi(i), wp) ! put in a complex number
        end do
        q = reverse(p)      !
        cp_ = reversez(cp)  ! for the ones that require reverse order

        write(*,'(A,1X,I3)')          ' Degree: ', degree
        write(*,'(A,1X/,*(g23.15/))') ' Coefficients: ', p(1:degree+1)

        if (icase==90 .or. icase==113) then
            write(*,*) 'skipping this case'
            cycle
        end if

        if (degree==2) then
            ! also test this one (only for quadratic equations):
            call dqdcrt(q, zr, zi)
            call check_results('dqdcrt', 0, zr, zi, degree)
        end if

        if (degree==3) then
            ! also test these (only for cubic equations):
            call dcbcrt(q, zr, zi)
            call check_results('dcbcrt', 0, zr, zi, degree)

            call rroots_chebyshev_cubic(p, zr, zi)
            call check_results('rroots_chebyshev_cubic', 0, zr, zi, degree)
        end if

        if (wp /= real128) then
            call polyroots(degree, p, zr, zi, istatus)
            call check_results('polyroots', istatus, zr, zi, degree)

            call cpolyroots(degree, cp, r, istatus)
            call check_results_complex('cpolyroots [complex coefficients]', istatus, real(r, wp), aimag(r), degree)
        end if

        call rpoly(p, degree, zr, zi, istatus)
        call check_results('rpoly', istatus, zr, zi, degree)

        istatus = 0 ! no estimates input
        call rpzero(degree,p,r,istatus,s)
        call check_results('rpzero', istatus, real(r,wp), aimag(r), degree)

        call rpqr79(degree,p,r,istatus)
        call check_results('rpqr79', istatus, real(r,wp), aimag(r), degree)

        call dpolz(degree,p,zr,zi,istatus)
        call check_results('dpolz', istatus, zr, zi, degree)

        call cpolz(cp,degree,r,istatus)
        call check_results_complex('cpolz [complex coefficients]', istatus, real(r,wp), aimag(r), degree)

        istatus = 0
        call cpoly(p,pi,degree,zr,zi,fail)
        if (fail) istatus = -1
        call check_results_complex('cpoly [complex coefficients]', istatus, zr, zi, degree)

        call cpqr79(degree,cp,r,istatus)
        call check_results_complex('cpqr79 [complex coefficients]', istatus, real(r,wp), aimag(r), degree)

        call qr_algeq_solver(degree,p,zr,zi,istatus,detil=detil)
        call check_results('qr_algeq_solver', istatus, zr,zi, degree)

        ! ,... or add to fpm.toml
        ![dev-dependencies]
        !    eiscor = { git="https://github.com/jacobwilliams/eiscor.git" }
        ! if (degree >=2) then
        !     call z_poly_roots(degree,cp,r,zr,istatus) ! just use zr for the residuals
        !     call check_results_complex('z_poly_roots [complex coefficients]', istatus, real(r, wp), aimag(r), degree)
        ! end if

        !....
        ! these accept the complex coefficients in reverse order
        call cmplx_roots_gen(degree, cp_, r) ! no status flag
        call check_results_complex('cmplx_roots_gen [complex coefficients]', 0, real(r, wp), aimag(r), degree)

        call polzeros(degree, cp_, 100, r, radius, err); istatus = 0; if (any(err)) istatus = -1
        call check_results_complex('polzeros [complex coefficients]', istatus, real(r, wp), aimag(r), degree)

        call fpml(cp_, degree, r, berr, cond, conv, itmax=100)
        call check_results_complex('fpml [complex coefficients]', 0, real(r, wp), aimag(r), degree)
        !....
        if (failure) error stop 'At least one test failed'
    end do

    !if (failure) error stop 'At least one test failed'

    contains

    !********************************************************************
        pure function reverse(x) result(y)

        !! reverse a `real(wp)` vector

        implicit none

        real(wp), dimension(:), intent(in) :: x
        real(wp), dimension(size(x)) :: y

        integer :: i !! counter
        integer :: n !! size of `x`

        n = size(x)

        do i = 1, n
            y(i) = x(n-i+1)
        end do

        end function reverse
    !********************************************************************

    !********************************************************************
        pure function reversez(x) result(y)

        !! reverse a `complex(wp)` vector

        implicit none

        complex(wp), dimension(:), intent(in) :: x
        complex(wp), dimension(size(x)) :: y

        integer :: i !! counter
        integer :: n !! size of `x`

        n = size(x)

        do i = 1, n
            y(i) = x(n-i+1)
        end do

        end function reversez
    !********************************************************************

    !********************************************************************
        subroutine allocate_arrays(d)

        integer,intent(in) :: d

        integer :: i

        degree = d

        p        = [(0, i=1,degree+1)]
        pi       = [(0, i=1,degree+1)]
        q        = [(0, i=1,degree+1)]
        cp       = [(0, i=1,degree+1)]
        cp_      = [(0, i=1,degree+1)]
        berr     = [(0, i=1,degree+1)]

        zr       = [(0, i=1,degree)]
        zi       = [(0, i=1,degree)]
        s        = [(0, i=1,degree)]
        r        = [(0, i=1,degree)]
        radius   = [(0, i=1,degree)]
        err      = [(.false., i=1,degree)]
        cond     = [(0, i=1,degree)]
        conv     = [(0, i=1,degree)]

        end subroutine allocate_arrays
    !********************************************************************

    !*****************************************************************************************
        subroutine check_results(name, istatus, zr, zi, degree)

            !! check the results.
            !! if any are not within the tolerance,
            !! then also try to polish them using the newton method.

            character(len=*),intent(in) :: name !! name of method
            integer,intent(in) :: istatus !! status flag (0 = success)
            real(wp),dimension(:),intent(in) :: zr, zi
            integer,intent(in) :: degree

            real(wp) :: zr_, zi_ ! copy of inputs for polishing
            real(wp),dimension(size(zr)) :: re, im ! copy of inputs for sorting
            complex(wp) :: z, root
            integer :: i,j !! counter
            integer :: istat

            real(wp),parameter :: tol = 1.0e-2_wp  !! acceptable root tolerance for tests
            real(wp),parameter :: ftol = 1.0e-8_wp !! desired root tolerance
            real(wp),parameter :: ztol = 10*epsilon(1.0_wp) !! newton tol for x
            logical,parameter :: polish = .true.

            write(*, '(/A,1x,i3)') trim(name)

            if (istatus /= 0) then
                failure = .true.
                write(*,'(A,1x,i3)') 'Error: method did not converge. istatus = ', istatus
                !error stop 'Error: method did not converge'
                return
            end if

            ! sort them in increasing order:
            re = zr
            im = zi
            call sort_roots(re, im)

            write(*, '(a)') '  real part               imaginary part         root'

            do j = 1, degree
                z = cmplx(re(j), im(j), wp)
                root = p(1)
                do i = 2, degree+1
                    root = root * z + p(i) ! horner's rule
                end do
                write(*, '(3(2g23.15,1x))') re(j), im(j), abs(root)
                if (polish .and. abs(root) > ftol) then
                    ! attempt to polish the root:
                    zr_ = re(j)
                    zi_ = im(j)
                    call newton_root_polish(degree, p, zr_, zi_, &
                                            ftol=ftol, ztol=ztol, maxiter=10, &
                                            istat=istat)
                    z = cmplx(zr_, zi_, wp) ! recompute root with possibly updated values
                    root = p(1)
                    do i = 2, degree+1
                        root = root * z + p(i) ! horner's rule
                    end do
                    write(*, '(3(2g23.15,1x),1X,A)') zr_, zi_, abs(root), 'POLISHED'
                    if (abs(root) > tol) then
                        failure = .true.
                        write(*,'(A)') 'Error: insufficient accuracy *******'
                        error stop 'Error: insufficient accuracy'
                    end if
                end if
            end do

        end subroutine check_results
    !*****************************************************************************************

    !*****************************************************************************************
        subroutine check_results_complex(name, istatus, zr, zi, degree)

            !! check the results (for complex coefficients).
            !! if any are not within the tolerance,
            !! then also try to polish them using the newton method.

            character(len=*),intent(in) :: name !! name of method
            integer,intent(in) :: istatus !! status flag (0 = success)
            real(wp),dimension(:),intent(in) :: zr, zi
            integer,intent(in) :: degree

            real(wp) :: zr_, zi_ ! copy of inputs for polishing
            real(wp),dimension(size(zr)) :: re, im ! copy of inputs for sorting
            complex(wp) :: z, root
            integer :: i,j !! counter
            integer :: istat

            real(wp),parameter :: tol = 1.0e-2_wp  !! acceptable root tolerance for tests
            real(wp),parameter :: ftol = 1.0e-8_wp !! desired root tolerance
            real(wp),parameter :: ztol = 10*epsilon(1.0_wp) !! newton tol for x
            logical,parameter :: polish = .true.

            write(*, '(/A,1x,i3)') trim(name)

            if (istatus /= 0) then
                failure = .true.
                write(*,'(A,1x,i3)') 'Error: method did not converge. istatus = ', istatus
                return
            end if

            ! sort them in increasing order:
            re = zr
            im = zi
            call sort_roots(re, im)

            write(*, '(a)') '  real part               imaginary part         root'

            do j = 1, degree
                z = cmplx(re(j), im(j), wp)
                root = cp(1)
                do i = 2, degree+1
                    root = root * z + cp(i) ! horner's rule
                end do
                write(*, '(3(2g23.15,1x))') re(j), im(j), abs(root)
                if (polish .and. abs(root) > ftol) then
                    ! attempt to polish the root:
                    zr_ = re(j)
                    zi_ = im(j)
                    call newton_root_polish(degree, cp, zr_, zi_, &
                                            ftol=ftol, ztol=ztol, maxiter=10, &
                                            istat=istat)
                    z = cmplx(zr_, zi_, wp) ! recompute root with possibly updated values
                    root = cp(1)
                    do i = 2, degree+1
                        root = root * z + cp(i) ! horner's rule
                    end do
                    write(*, '(3(2g23.15,1x),1X,A)') zr_, zi_, abs(root), 'POLISHED'
                    if (abs(root) > tol) then
                        failure = .true.
                        write(*,'(A)') 'Error: insufficient accuracy *******'
                        error stop 'Error: insufficient accuracy'
                    end if
                end if
            end do

        end subroutine check_results_complex
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

            !call random_number(x)
            x = rand%genrand64_real1()

            x = a + (b-a)*x

        end function get_random_number
    !*****************************************************************************************

!*****************************************************************************************
    end program polyroots_test
!*****************************************************************************************