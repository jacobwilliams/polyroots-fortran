!*****************************************************************************************
!>
!  Compute roots for all 10th degree polynomials with coefficients +/- 1.

    program polyroots_test_10

    use polyroots_module, only: dpolz, wp => polyroots_module_rk
    use pyplot_module,    only: pyplot

    implicit none

    integer,parameter :: degree = 10 !! polynomial degree
    integer,dimension(*),parameter :: icoeffs = [-1,1] !! set of coefficients
    integer,dimension (degree+1) :: a !! coefficients of polynomial
    real(wp),dimension(degree) :: zr, zi !! roots
    integer :: ierr !! error code from [[dpolz]]
    type(pyplot) :: plt !! for making the plot

    call plt%initialize(grid=.true.,xlabel='$\Re(z)$',ylabel='$\Im(z)$',&
                        title='Degree 10 Polynomial Roots',usetex=.true.,&
                        figsize=[20,10])

    call generate(1)
    call plt%savefig('roots.png')

    contains

        recursive subroutine generate (i)
        integer, intent(in) :: i
        integer :: ix
        if (i > degree+1) then
            !write (*, '(*(I2,","))') a
            call dpolz(degree,real(a,wp),zr,zi,ierr); if (ierr/=0) error stop ierr
            call plt%add_plot(zr,zi,label='',linestyle='bo',markersize=1)
        else
            do ix = 1,size(icoeffs) 
                a(i) = icoeffs(ix)
                call generate(i+1)
            end do
        end if
        end subroutine generate

    end program polyroots_test_10
!*****************************************************************************************