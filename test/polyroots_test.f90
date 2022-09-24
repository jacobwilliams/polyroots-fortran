!*****************************************************************************************
!>
!  Tests for [[rpoly]].

    program polyroots_test

    use iso_fortran_env
    use polyroots_module, wp => polyroots_module_rk

    implicit none

    integer,parameter :: max_degree = 10 !! max degree polynomials to test for random cases
    integer,parameter :: n_cases = 14 + 110 !! number of cases to run

    real(wp),dimension(:),allocatable :: p, zr, zi, s, q, radius,rr,rc, berr,cond
    integer,dimension(:),allocatable :: conv
    complex(wp),dimension(:),allocatable :: r, cp
    integer :: degree, i, istatus, icase, n
    integer,dimension(:),allocatable :: seed
    real(wp) :: detil
    logical :: fail
    logical :: failure !! if any of the tests failed
    logical,dimension(:),allocatable :: err
    real(wp) :: x1,x2,x3
    integer :: l
    integer :: idegree !! counter for degrees to test
    integer :: n_degree !! number of tests run for each degree so far

    failure = .false.

    ! set random seed for consistent results:
    call random_seed(size=n)
    allocate(seed(n))
    seed = 42
    call random_seed(put=seed)
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
        case(2)
            call allocate_arrays(4)
            p = [1,-3,20,44,54]
        case(3)
            call allocate_arrays(6)
            p = [1,-2,2,1,6,-6,8]
        case(4)
            call allocate_arrays(5)
            p = [1,1,-8,-16,7,15]
        case(5)
            call allocate_arrays(5)
            p = [1,7,5,6,3,2]
        case(6)
            call allocate_arrays(5)
            p = [2,3,6,5,7,1]
        case(7)
            call allocate_arrays(6)
            p = [1,0,-14,0,49,0,-36]
        case(8)
            call allocate_arrays(8)
            p = [1,0,-30,0,273,0,-820,0,576]
        case(9)
            call allocate_arrays(4)
            p = [1,0,0,0,-16]
        case(10)
            call allocate_arrays(6)
            p = [1,-2,2,1,6,-6,8]
        case(11)
            ! a case where 1 is an obvious root
            call allocate_arrays(5)
            p = [8,-8,16,-16,8,-8]
        case(12)
            call allocate_arrays(3)
            p = [ -8.0e18_wp,3.0e12_wp,5.0e6_wp,1.0_wp]
        case(13)
            call allocate_arrays(3)
            p = [4.0_wp, 3.0_wp, 2.0_wp, 1.0_wp]
        case(14)
            call allocate_arrays(2)
            p = [3.0_wp, 2.0_wp, 1.0_wp]

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
            end do
        end select

        write(*,'(A,1X,I3)')          ' Degree: ', degree
        write(*,'(A,1X/,*(g23.15/))') ' Coefficients: ', p(1:degree+1)

        q = reverse(p) ! the following two accept the coefficients in reverse order
        do i = 1, degree+1
            cp(i) = cmplx(p(i), 0.0_wp, wp) ! put in a complex number
        end do

        if (degree==2) then
            ! also test this one (only for quadratic equations):
            call dqdcrt(q, zr, zi)
            call check_results('dqdcrt', 0, zr, zi, degree)
        end if

        if (degree==3) then
            ! also test these (only for cubic equations):
            call dcbcrt(q, zr, zi)
            call check_results('dcbcrt', 0, zr, zi, degree)

            call lebedev(p, zr, zi)
            call check_results('lebedev', 0, zr, zi, degree)
        end if

        if (wp /= REAL128) then
            call polyroots(degree, p, zr, zi, istatus)
            call check_results('polyroots', istatus, zr, zi, degree)

            call cpolyroots(degree, cp, r, istatus)
            call check_results('cpolyroots', istatus, real(r, wp), aimag(r), degree)
        end if

        call rpoly(p, degree, zr, zi, istatus)
        call check_results('rpoly', istatus, zr, zi, degree)

        istatus = 0 ! no estimates input
        call rpzero(degree,p,r,istatus,s)
        call check_results('rpzero', istatus, real(r,wp), aimag(r), degree)

        call rpqr79(degree,p,r,istatus)
        call check_results('rpqr79', istatus, real(r,wp), aimag(r), degree)

        ! for now, just test the following two with the real coefficients only:

        q = 0.0_wp
        istatus = 0
        call cpoly(p,q,degree,zr,zi,fail)
        if (fail) istatus = -1
        call check_results('cpoly', istatus, zr, zi, degree)

        call cpqr79(degree,cp,r,istatus)
        call check_results('cpqr79', istatus, real(r,wp), aimag(r), degree)

        call qr_algeq_solver(degree,p,zr,zi,istatus,detil=detil)
        call check_results('qr_algeq_solver', istatus, real(r,wp), aimag(r), degree)

        !..... these accept the complex coefficients in reverse order
        cp = reversez(cp)
        call cmplx_roots_gen(degree, cp, r) ! no status flag
        rr = real(r, wp)
        rc = aimag(r)
        call check_results('cmplx_roots_gen', 0, rr, rc, degree)

        istatus = 0
        call polzeros(degree, cp, 100, r, radius, err)
        if (any(err)) istatus = -1
        rr = real(r, wp)
        rc = aimag(r)
        call check_results('polzeros', istatus, rr, rc, degree)

        call fpml(cp, degree, r, berr, cond, conv, itmax=100)
        call check_results('fpml', 0, real(r, wp), aimag(r), degree)

    end do

    if (failure) error stop 'At least one test failed'

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
        q        = [(0, i=1,degree+1)]
        cp       = [(0, i=1,degree+1)]
        berr     = [(0, i=1,degree+1)]

        zr       = [(0, i=1,degree)]
        zi       = [(0, i=1,degree)]
        s        = [(0, i=1,degree)]
        r        = [(0, i=1,degree)]
        radius   = [(0, i=1,degree)]
        err      = [(.false., i=1,degree)]
        rr       = [(0, i=1,degree)]
        rc       = [(0, i=1,degree)]
        cond     = [(0, i=1,degree)]
        conv     = [(0, i=1,degree)]

        end subroutine allocate_arrays
    !********************************************************************

    !*****************************************************************************************
        subroutine check_results(name, istatus, re, im, degree)

            !! check the results.
            !! if any are not within the tolerance,
            !! then also try to polish them using the newton method.

            character(len=*),intent(in) :: name !! name of method
            integer,intent(in) :: istatus !! status flag (0 = success)
            real(wp),dimension(:),intent(in) :: re, im
            integer,intent(in) :: degree

            complex(wp) :: z, root
            integer :: i,j !! counter
            real(wp) :: zr, zi ! copy of inputs for polishing
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
                    if (abs(root) > tol) then
                        failure = .true.
                        write(*,'(A)') 'Error: insufficient accuracy *******'
                        !error stop 'Error: insufficient accuracy'
                    end if
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
    end program polyroots_test
!*****************************************************************************************