!*****************************************************************************************
!>
!  Tests for [[rpoly]].

    program polyroots_test

    use iso_fortran_env
    use polyroots_module, wp => polyroots_module_rk

    implicit none

    integer,parameter :: max_degree = 10 !! max degree polynomials to test for random cases
    integer,parameter :: n_cases = 30 !! number of cases to run

    real(wp),dimension(:),allocatable :: p, zr, zi, s, q, radius
    complex(wp),dimension(:),allocatable :: r, cp
    integer :: degree, i, istatus, icase, n
    integer,dimension(:),allocatable :: seed
    real(wp) :: detil
    logical :: fail
    logical :: failure !! if any of the tests failed
    logical,dimension(:),allocatable :: err

    !--------------------------------------

    failure = .false.

    ! set random seed for consistent results:
    call random_seed(size=n)
    allocate(seed(n))
    seed = 42
    call random_seed(put=seed)

    do icase = 1, 1 !n_cases

        write(*,'(/A,I2,A)') '--------CASE ', icase, ' ---------'

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
            ! random coefficients
            call allocate_arrays(get_random_integer_number(3,max_degree))
            do i = 1, degree+1
                p(i)  = get_random_number(-1000.0_wp,1000.0_wp)
            end do
        end select

        write(*,'(A,1X,I3)')          ' Degree: ', degree
        write(*,'(A,1X/,*(g23.15/))') ' Coefficients: ', p(1:degree+1)

        q = reverse(p) ! the following two accept the coefficients in reverse order

        ! if (degree==2) then
        !     ! also test this one (only for quadratic equations):
        !     write(*, '(A,1x,i3)') 'dqdcrt'
        !     write(*, '(a)') '  real part               imaginary part         root'
        !     call dqdcrt(q, zr, zi)
        !     call check_results(0, zr, zi, degree)
        ! end if

        ! if (degree==3) then
        !     ! also test this one (only for cubic equations):
        !     write(*, '(A,1x,i3)') 'dcbcrt'
        !     write(*, '(a)') '  real part               imaginary part         root'
        !     call dcbcrt(q, zr, zi)
        !     call check_results(0, zr, zi, degree)
        ! end if

        ! write(*, '(A,1x,i3)') 'rpoly'
        ! write(*, '(a)') '  real part               imaginary part         root'
        ! call rpoly(p, degree, zr, zi, istatus)
        ! call check_results(istatus, zr, zi, degree)

        ! write(*, '(/A,1x,i3)') 'rpzero'
        ! write(*, '(a)') '  real part               imaginary part         root'
        ! istatus = 0 ! no estimates input
        ! call rpzero(degree,p,r,istatus,s)
        ! call check_results(istatus,real(r,wp), aimag(r), degree)

        ! write(*, '(/A,1x,i3)') 'rpqr79'
        ! write(*, '(a)') '  real part               imaginary part         root'
        ! call rpqr79(degree,p,r,istatus)
        ! call check_results(istatus,real(r,wp), aimag(r), degree)

        ! ! for now, just test the following two with the real coefficients only:

        ! write(*, '(/A,1x,i3)') 'cpoly'
        ! write(*, '(a)') '  real part               imaginary part         root'
        ! q = 0.0_wp
        ! istatus = 0
        ! call cpoly(p,q,degree,zr,zi,fail)
        ! if (fail) istatus = -1
        ! call check_results(istatus, zr, zi, degree)

        ! write(*, '(/A,1x,i3)') 'cpqr79'
        ! write(*, '(a)') '  real part               imaginary part         root'
        do i = 1, degree+1
            cp(i) = cmplx(p(i), 0.0_wp, wp) ! put in a complex number
        end do
        ! call cpqr79(degree,cp,r,istatus)
        ! call check_results(istatus, real(r,wp), aimag(r), degree)

        ! write(*, '(/A,1x,i3)') 'qr_algeq_solver'
        ! write(*, '(a)') '  real part               imaginary part         root'
        ! call qr_algeq_solver(degree,p,zr,zi,istatus,detil=detil)
        ! call check_results(istatus, real(r,wp), aimag(r), degree)

        ! !..... these accept the complex coefficients in reverse order
        ! write(*, '(/A,1x,i3)') 'cmplx_roots_gen'
        ! write(*, '(a)') '  real part               imaginary part         root'
        cp = reversez(cp)
        ! call cmplx_roots_gen(degree, cp, r) ! no status flag
        ! call check_results(0, zr, zi, degree)

        write(*, '(/A,1x,i3)') 'polzeros'
        write(*, '(a)') '  real part               imaginary part         root'
        istatus = 0
        call polzeros(degree, cp, 100, r, radius, err)
        if (any(err)) istatus = -1
        call check_results(istatus, real(r, wp), aimag(r), degree)

        ! if (wp /= REAL128) then
        !     write(*, '(/A,1x,i3)') 'polyroots'
        !     write(*, '(a)') '  real part               imaginary part         root'
        !     call polyroots(degree, p, zr, zi, istatus)
        !     call check_results(istatus, zr, zi, degree)
        ! end if

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

        degree = d

        if (allocated(p))  deallocate(p)
        if (allocated(zr)) deallocate(zr)
        if (allocated(zi)) deallocate(zi)
        if (allocated(s))  deallocate(s)
        if (allocated(r))  deallocate(r)
        if (allocated(cp)) deallocate(cp)
        if (allocated(q))  deallocate(q)
        if (allocated(radius))  deallocate(radius)
        if (allocated(err))  deallocate(err)

        allocate(p(degree+1))
        allocate(q(degree+1))
        allocate(zr(degree+1))
        allocate(zi(degree+1))
        allocate(s(degree+1))
        allocate(r(degree+1))
        allocate(cp(degree+1))
        allocate(radius(degree))
        allocate(err(degree))

        end subroutine allocate_arrays
    !********************************************************************

    !*****************************************************************************************
        subroutine check_results(istatus, re, im, degree)

            !! check the results.
            !! if any are not within the tolerance,
            !! then also try to polish them using the newton method.

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
            logical,parameter :: polish = .false. !! skip polishing

            if (istatus /= 0) then
                failure = .true.
                write(*,'(A,1x,i3)') 'Error: method did not converge. istatus = ', istatus
                return
            end if

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

    subroutine sort(x)

    !! sort the vector x, according to nonincreasing real parts

    complex(wp),dimension(:),intent(inout) :: x

    integer     :: n,k,i,imax
    complex(wp) :: temp
    real(wp)    :: amax,rxi

    n = size(x)

    do k=1,n-1
        amax=real(x(k), wp)
        imax=k
        do i=k+1,n
            rxi=real(x(i), wp)
            if (amax<rxi) then
                amax=rxi
                imax=i
            endif
        end do
        temp=x(k)
        x(k)=x(imax)
        x(imax)=temp
    end do

    end subroutine sort

!*****************************************************************************************
    end program polyroots_test
!*****************************************************************************************