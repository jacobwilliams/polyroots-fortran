!*****************************************************************************************
!>
!  Example in the readme.

program example

use iso_fortran_env
use polyroots_module, wp => polyroots_module_rk

implicit none

integer,parameter :: degree = 5 !! polynomial degree
real(wp),dimension(degree+1) :: p = [1,2,3,4,5,6] !! coefficients

integer :: i !! counter
integer :: istatus !! status code
real(wp),dimension(degree) :: zr !! real components of roots
real(wp),dimension(degree) :: zi !! imaginary components of roots

call polyroots(degree, p, zr, zi, istatus)

write(*,'(A,1x,I3)') 'istatus: ', istatus
write(*, '(*(a22,1x))') 'real part', 'imaginary part'
do i = 1, degree
    write(*,'(*(e22.15,1x))') zr(i), zi(i)
end do

end program example
!*****************************************************************************************