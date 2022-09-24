!*****************************************************************************************
!>
!  Tests for [[dcbcrt]] and [[dqdcrt]].

    program dcbcrt_test

    use iso_fortran_env
    use polyroots_module, wp => polyroots_module_rk

    implicit none

    real(wp), dimension(4)  :: a    !! coefficients
    real(wp), dimension(3)  :: zr   !! real components of roots
    real(wp), dimension(3)  :: zi   !! imaginary components of roots

    integer :: i !! counter
    complex(wp) :: z, root

    a = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp]

    write(*,'(/A)') 'dcbcrt test:'
    call dcbcrt(a, zr, zi)
    do i = 1, 3
        z = cmplx(zr(i), zi(i), wp)
        root = a(1) + a(2)*z + a(3)*z**2 + a(4)*z**3
        write(*,'(A,1x,*(e22.15,1x))') 'root is: ', zr(i), zi(i), abs(root)
        if (abs(root) > 100*epsilon(1.0_wp)) error stop 'Error: insufficient accuracy'
    end do

    write(*,'(/A)') 'dqdcrt test:'
    call dqdcrt(a(1:3), zr(1:2), zi(1:2))
    do i = 1, 2
        z = cmplx(zr(i), zi(i), wp)
        root = a(1) + a(2)*z + a(3)*z**2
        write(*,'(A,1x,*(e22.15,1x))') 'root is: ', zr(i), zi(i), abs(root)
        if (abs(root) > 100*epsilon(1.0_wp)) error stop 'Error: insufficient accuracy'
    end do

    end program dcbcrt_test
!*****************************************************************************************