!*****************************************************************************************
!>
!  Compute roots for all 11th degree polynomials with coefficients +/- 1.

    program polyroots_test_11

    use polyroots_module, only: polyroots, wp => polyroots_module_rk
    use pyplot_module,    only: pyplot

    implicit none

    integer,parameter :: degree = 11 !! polynomial degree

    integer,dimension(2) :: icoeffs = [-1,1] !! set of coefficients
    integer :: ierr !! error code from [[polyroots]]
    type(pyplot) :: plt !! for making the plot
    integer :: j !! counter
    integer,dimension(degree+1) :: a !! coefficients of polynomial

    call plt%initialize(grid=.true.,xlabel='$\\Re(z)$',ylabel='$\\Im(z)$',&
                        title='Degree 11 Polynomial Roots',usetex=.true.,&
                        font_size=30, axes_labelsize=30, &
                        xtick_labelsize=30, ytick_labelsize=30, &
                        figsize=[20,10])
    call generate(1)
    call plt%savefig('roots_11.png')
    ! call plt%showfig()

    contains

        recursive subroutine generate (i)
        integer, intent(in) :: i
        integer :: ix
        integer,save :: iprint = 0
        integer :: k
        real(wp),dimension(degree) :: zr, zi !! roots
        if (i > degree+1) then
            iprint = iprint + 1
            if (iprint==1) write (*, '(*(I2,","))') a
            call polyroots(degree,real(a,wp),zr,zi,ierr)
            if (iprint==1) then
                ! print the first set of roots
                do k = 1, degree
                    write(*,'(ES20.12, a, ES20.12, a)') zr(k), ' + ', zi(k), ' i'
                    ! write(*,'(A,I2,A,ES12.5)') 'zr(',k,') = ', zr(k)
                    ! write(*,'(A,I2,A,ES12.5)') 'zi(',k,') = ', zi(k)
                end do
            end if
            if (ierr/=0) error stop 'error calling polyroots'
            call plt%add_plot(zr,zi,label='',linestyle='bo',markersize=2)
        else
            do ix = 1,size(icoeffs)
                a(i) = icoeffs(ix)
                call generate(i+1)
            end do
        end if
        end subroutine generate

    end program polyroots_test_11
!*****************************************************************************************