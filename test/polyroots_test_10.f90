!*****************************************************************************************
!>
!  Compute roots for all 10th degree polynomials with coefficients +/- 1.

    program polyroots_test_10

    use polyroots_module, only: polyroots, dpolz, wp => polyroots_module_rk
    use pyplot_module,    only: pyplot

    implicit none

    integer :: degree = 5 !! polynomial degree
    integer,dimension(2) :: icoeffs = [-1,1] !! set of coefficients
    integer :: ierr !! error code from [[dpolz]]
    type(pyplot) :: plt !! for making the plot
    integer :: i ,j
    character(len=5) :: istr
    integer,dimension (:),allocatable :: a !! coefficients of polynomial
    real(wp),dimension(:),allocatable :: zr, zi !! roots

    do i = 2, 10, 2

        write(istr,'(I5)') i; istr = adjustl(istr)
        degree = i
        ! resize
        a  = [(0, j = 1, degree+1)]
        zr = [(0, j = 1, degree)]
        zi = [(0, j = 1, degree)]

        !icoeffs = icoeffs + 1
        write(*,*) 'degree = ', degree
        call plt%initialize(grid=.true.,xlabel='$\Re(z)$',ylabel='$\Im(z)$',&
                            title='Degree '//trim(istr)//' Polynomial Roots',usetex=.true.,&
                            font_size=25, axes_labelsize=25, &
                            xtick_labelsize=25, ytick_labelsize=25, &
                            figsize=[20,10])
        call generate(1)
        call plt%savefig('roots_'//trim(istr)//'.png')

    end do

    contains

        recursive subroutine generate (i)
        integer, intent(in) :: i
        integer :: ix
        if (i > degree+1) then
            !write (*, '(*(I2,","))') a
            call polyroots(degree,real(a,wp),zr,zi,ierr)  !polyroots !! dpolz
            if (ierr/=0) return !error stop ierr
            call plt%add_plot(zr,zi,label='',linestyle='bo',markersize=2)
        else
            do ix = 1,size(icoeffs)
                a(i) = icoeffs(ix)
                call generate(i+1)
            end do
        end if
        end subroutine generate

    end program polyroots_test_10
!*****************************************************************************************