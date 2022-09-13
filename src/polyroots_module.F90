!*****************************************************************************************
!>
!  Polynomial Roots with Modern Fortran
!
!### Author
!  * Jacob Williams
!
!@note The default real kind (`wp`) can be
!      changed using optional preprocessor flags.
!      This library was built with real kind:
#ifdef REAL32
!      `real(kind=real32)` [4 bytes]
#elif REAL64
!      `real(kind=real64)` [8 bytes]
#elif REAL128
!      `real(kind=real128)` [16 bytes]
#else
!      `real(kind=real64)` [8 bytes]
#endif

module polyroots_module

   use iso_fortran_env

   implicit none

   private

#ifdef REAL32
   integer, parameter, public :: polyroots_module_rk = real32   !! real kind used by this module [4 bytes]
#elif REAL64
   integer, parameter, public :: polyroots_module_rk = real64   !! real kind used by this module [8 bytes]
#elif REAL128
   integer, parameter, public :: polyroots_module_rk = real128  !! real kind used by this module [16 bytes]
#else
   integer, parameter, public :: polyroots_module_rk = real64   !! real kind used by this module [8 bytes]
#endif

   integer, parameter :: wp = polyroots_module_rk  !! local copy of `polyroots_module_rk` with a shorter name

   real(wp), parameter :: eps = epsilon(1.0_wp)

   ! general polynomial routines:
   public :: rpoly

   ! special polynomial routines:
   public :: dcbcrt
   public :: dqdcrt

   ! utilities routines:
   public :: dcbrt

contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finds the zeros of a general real polynomial using the jenkins & traub algorithm
!
!### History
! * algorithm 493 by jenkins & traub
! * code converted using to_f90 by alan miller, 2003-06-02
! * Jacob Williams, 9/13/2022 : modernized this code

    subroutine rpoly(op, degree, zeror, zeroi, fail)

      implicit none

      real(wp), intent(in)    :: op(:) !! vector of coefficients in order of decreasing powers
      integer, intent(inout)  :: degree !! degree of polynomial
      real(wp), intent(out)   :: zeror(:) !! output vector of real parts of the zeros
      real(wp), intent(out)   :: zeroi(:) !! output vector of imaginary parts of the zeros
      logical, intent(out)    :: fail !! true only if leading coefficient is zero
                                      !! or if rpoly has found fewer than degree zeros.
                                      !! in the latter case degree is reset to the number of zeros found.

      ! these were formerly in a common block:
      real(wp), allocatable   :: p(:), qp(:), k(:), qk(:), svk(:)
      real(wp) :: sr, si, u, v, a, b, c, d, a1, a3, &
                                 a7, e, f, g, h, szr, szi, lzr, lzi
      integer :: n, nn
      !--------------------------------------------------------

      real(wp), allocatable :: temp(:)
      real(wp), allocatable :: pt(:)
      real(wp) :: t, aa, bb, cc, factor, mx, mn, xx, yy, xxx, x, sc, bnd, xm, ff, df, dx
      integer   :: cnt, nz, i, j, jj, l, nm1
      logical   :: zerok

      real(wp),parameter :: deg2rad = acos(-1.0_wp) / 180.0_wp
      real(wp),parameter :: cosr = cos(94.0_wp * deg2rad)
      real(wp),parameter :: sinr = sin(86.0_wp * deg2rad)
      real(wp),parameter :: base = radix(0.0_wp)
      real(wp),parameter :: eta = epsilon(1.0_wp)
      real(wp),parameter :: infin = huge(0.0_wp)
      real(wp),parameter :: smalno = tiny(0.0_wp)
      real(wp),parameter :: sqrthalf = sqrt(0.5_wp)
      real(wp),parameter :: are = eta !! unit error in +
      real(wp),parameter :: mre = eta !! unit error in *
      real(wp),parameter :: lo = smalno/eta

      ! initialization of constants for shift rotation
      xx = sqrthalf
      yy = -xx
      fail = .false.
      n = degree
      nn = n + 1

      ! algorithm fails if the leading coefficient is zero.
      if (op(1) == 0.0_wp) then
         fail = .true.
         degree = 0
         return
      end if

      ! remove the zeros at the origin if any
      do
         if (op(nn) /= 0.0_wp) exit
         j = degree - n + 1
         zeror(j) = 0.0_wp
         zeroi(j) = 0.0_wp
         nn = nn - 1
         n = n - 1
     end do

      ! allocate various arrays
      if (allocated(p)) deallocate (p, qp, k, qk, svk)
      allocate (p(nn), qp(nn), k(nn), qk(nn), svk(nn), temp(nn), pt(nn))

      ! make a copy of the coefficients
      p(1:nn) = op(1:nn)

      main : do

        ! start the algorithm for one zero
        if (n <= 2) then
            if (n < 1) return

            ! calculate the final zero or pair of zeros
            if (n /= 2) then
                zeror(degree) = -p(2)/p(1)
                zeroi(degree) = 0.0_wp
                return
            end if
            call quad(p(1), p(2), p(3), zeror(degree - 1), zeroi(degree - 1), &
                      zeror(degree), zeroi(degree))
            return
        end if

        ! find largest and smallest moduli of coefficients.
        mx = 0.0_wp ! max
        mn = infin  ! min
        do i = 1, nn
            x = abs(real(p(i), wp))
            if (x > mx) mx = x
            if (x /= 0.0_wp .and. x < mn) mn = x
        end do

        ! scale if there are large or very small coefficients computes a scale
        ! factor to multiply the coefficients of the polynomial.
        ! the scaling is done to avoid overflow and to avoid undetected underflow
        ! interfering with the convergence criterion.
        ! the factor is a power of the base
        scale : block
            sc = lo/mn
            if (sc <= 1.0_wp) then
                if (mx < 10.0_wp) exit scale
                if (sc == 0.0_wp) sc = smalno
            else
                if (infin/sc < mx) exit scale
            end if
            l = log(sc)/log(base) + 0.5_wp
            factor = (base*1.0_wp)**l
            if (factor /= 1.0_wp) then
                p(1:nn) = factor*p(1:nn)
            end if
        end block scale

        ! compute lower bound on moduli of zeros.
        pt(1:nn) = abs(p(1:nn))
        pt(nn) = -pt(nn)

        ! compute upper estimate of bound
        x = exp((log(-pt(nn)) - log(pt(1)))/n)
        if (pt(n) /= 0.0_wp) then
            ! if newton step at the origin is better, use it.
            xm = -pt(nn)/pt(n)
            if (xm < x) x = xm
        end if

        ! chop the interval (0,x) until ff <= 0
        do
            xm = x*0.1_wp
            ff = pt(1)
            do i = 2, nn
                ff = ff*xm + pt(i)
            end do
            if (ff > 0.0_wp) then
                x = xm
            else
                exit
            end if
        end do
        dx = x

        ! do newton iteration until x converges to two decimal places
        do
            if (abs(dx/x) <= 0.005_wp) exit
            ff = pt(1)
            df = ff
            do i = 2, n
                ff = ff*x + pt(i)
                df = df*x + ff
            end do
            ff = ff*x + pt(nn)
            dx = ff/df
            x = x - dx
        end do
        bnd = x

        ! compute the derivative as the intial k polynomial
        ! and do 5 steps with no shift
        nm1 = n - 1
        do i = 2, n
            k(i) = (nn - i)*p(i)/n
        end do
        k(1) = p(1)
        aa = p(nn)
        bb = p(n)
        zerok = k(n) == 0.0_wp
        do jj = 1, 5
            cc = k(n)
            if (.not. zerok) then
                ! use scaled form of recurrence if value of k at 0 is nonzero
                t = -aa/cc
                do i = 1, nm1
                j = nn - i
                k(j) = t*k(j - 1) + p(j)
                end do
                k(1) = p(1)
                zerok = abs(k(n)) <= abs(bb)*eta*10.0_wp
            else
                ! use unscaled form of recurrence
                do i = 1, nm1
                j = nn - i
                k(j) = k(j - 1)
                end do
                k(1) = 0.0_wp
                zerok = k(n) == 0.0_wp
            end if
        end do

        ! save k for restarts with new shifts
        temp(1:n) = k(1:n)

        ! loop to select the quadratic  corresponding to each
        ! new shift
        do cnt = 1, 20
            ! quadratic corresponds to a double shift to a non-real point and its complex
            ! conjugate.  the point has modulus bnd and amplitude rotated by 94 degrees
            ! from the previous shift
            xxx = cosr*xx - sinr*yy
            yy = sinr*xx + cosr*yy
            xx = xxx
            sr = bnd*xx
            si = bnd*yy
            u = -2.0_wp*sr
            v = bnd

            ! second stage calculation, fixed quadratic
            call fxshfr(20*cnt, nz)
            if (nz /= 0) then
                ! the second stage jumps directly to one of the third stage iterations and
                ! returns here if successful.
                ! deflate the polynomial, store the zero or zeros and return to the main
                ! algorithm.
                j = degree - n + 1
                zeror(j) = szr
                zeroi(j) = szi
                nn = nn - nz
                n = nn - 1
                p(1:nn) = qp(1:nn)
                if (nz /= 1) then
                    zeror(j + 1) = lzr
                    zeroi(j + 1) = lzi
                end if
                cycle main
            end if

            ! if the iteration is unsuccessful another quadratic
            ! is chosen after restoring k
            k(1:nn) = temp(1:nn)
        end do

        exit
      end do main

      ! return with failure if no convergence with 20 shifts
      fail = .true.
      degree = degree - n

   contains

      subroutine fxshfr(l2, nz)

         !! computes up to  l2  fixed shift k-polynomials, testing for convergence in
         !! the linear or quadratic case.  initiates one of the variable shift
         !! iterations and returns with the number of zeros found.

         integer, intent(in)   :: l2 !! limit of fixed shift steps
         integer, intent(out)  :: nz !! number of zeros found

         real(wp) :: svu, svv, ui, vi, s, betas, betav, oss, ovv, &
                     ss, vv, ts, tv, ots, otv, tvv, tss
         integer  :: type, j, iflag
         logical  :: vpass, spass, vtry, stry

         nz = 0
         betav = 0.25_wp
         betas = 0.25_wp
         oss = sr
         ovv = v

         ! evaluate polynomial by synthetic division
         call quadsd(nn, u, v, p, qp, a, b)
         call calcsc(type)
         do j = 1, l2
            ! calculate next k polynomial and estimate v
            call nextk(type)
            call calcsc(type)
            call newest(type, ui, vi)
            vv = vi

            ! estimate s
            ss = 0.0_wp
            if (k(n) /= 0.0_wp) ss = -p(nn)/k(n)
            tv = 1.0_wp
            ts = 1.0_wp
            if (j /= 1 .and. type /= 3) then
               ! compute relative measures of convergence of s and v sequences
               if (vv /= 0.0_wp) tv = abs((vv - ovv)/vv)
               if (ss /= 0.0_wp) ts = abs((ss - oss)/ss)

               ! if decreasing, multiply two most recent convergence measures
               tvv = 1.0_wp
               if (tv < otv) tvv = tv*otv
               tss = 1.0_wp
               if (ts < ots) tss = ts*ots

               ! compare with convergence criteria
               vpass = tvv < betav
               spass = tss < betas
               if (spass .or. vpass) then

                  ! at least one sequence has passed the convergence test.
                  ! store variables before iterating
                  svu = u
                  svv = v
                  svk(1:n) = k(1:n)
                  s = ss

                  ! choose iteration according to the fastest converging sequence
                  vtry = .false.
                  stry = .false.
                  if (spass .and. ((.not. vpass) .or. tss < tvv)) go to 40

20                call quadit(ui, vi, nz)
                  if (nz > 0) return

                  ! quadratic iteration has failed. flag that it has
                  ! been tried and decrease the convergence criterion.
                  vtry = .true.
                  betav = betav * 0.25_wp

                  ! try linear iteration if it has not been tried and
                  ! the s sequence is converging
                  if (stry .or. (.not. spass)) go to 50
                  k(1:n) = svk(1:n)
40                call realit(s, nz, iflag)
                  if (nz > 0) return

                  ! linear iteration has failed.  flag that it has been
                  ! tried and decrease the convergence criterion
                  stry = .true.
                  betas = betas * 0.25_wp
                  if (iflag /= 0) then
                     ! if linear iteration signals an almost double real
                     ! zero attempt quadratic interation
                     ui = -(s + s)
                     vi = s*s
                     go to 20
                  end if

                  ! restore variables
50                u = svu
                  v = svv
                  k(1:n) = svk(1:n)

                  ! try quadratic iteration if it has not been tried
                  ! and the v sequence is converging
                  if (vpass .and. (.not. vtry)) go to 20

                  ! recompute qp and scalar values to continue the second stage
                  call quadsd(nn, u, v, p, qp, a, b)
                  call calcsc(type)
               end if
            end if
            ovv = vv
            oss = ss
            otv = tv
            ots = ts
         end do

      end subroutine fxshfr

      subroutine quadit(uu, vv, nz)

         !! variable-shift k-polynomial iteration for a quadratic factor, converges
         !! only if the zeros are equimodular or nearly so.

         real(wp), intent(in)  :: uu !! coefficients of starting quadratic
         real(wp), intent(in)  :: vv !! coefficients of starting quadratic
         integer, intent(out)  :: nz !! number of zero found

         real(wp) :: ui, vi, mp, omp, ee, relstp, t, zm
         integer  :: type, i, j
         logical  :: tried

         nz = 0
         tried = .false.
         u = uu
         v = vv
         j = 0

         ! main loop
         main : do
            call quad(1.0_wp, u, v, szr, szi, lzr, lzi)

            ! return if roots of the quadratic are real and not
            ! close to multiple or nearly equal and  of opposite sign.
            if (abs(abs(szr) - abs(lzr)) > 0.01_wp*abs(lzr)) return

            ! evaluate polynomial by quadratic synthetic division
            call quadsd(nn, u, v, p, qp, a, b)
            mp = abs(a - szr*b) + abs(szi*b)

            ! compute a rigorous  bound on the rounding error in evaluting p
            zm = sqrt(abs(v))
            ee = 2.0_wp*abs(qp(1))
            t = -szr*b
            do i = 2, n
                ee = ee*zm + abs(qp(i))
            end do
            ee = ee*zm + abs(a + t)
            ee = (5.0_wp*mre + 4.0_wp*are)*ee - &
                 (5.0_wp*mre + 2.0_wp*are)*(abs(a + t) + &
                  abs(b)*zm) + 2.0_wp*are*abs(t)

            ! iteration has converged sufficiently if the
            ! polynomial value is less than 20 times this bound
            if (mp <= 20.0_wp*ee) then
                nz = 2
                return
            end if
            j = j + 1

            ! stop iteration after 20 steps
            if (j > 20) return
            if (j >= 2) then
                if (.not. (relstp > 0.01_wp .or. mp < omp .or. tried)) then

                ! a cluster appears to be stalling the convergence.
                ! five fixed shift steps are taken with a u,v close to the cluster
                if (relstp < eta) relstp = eta
                relstp = sqrt(relstp)
                u = u - u*relstp
                v = v + v*relstp
                call quadsd(nn, u, v, p, qp, a, b)
                do i = 1, 5
                    call calcsc(type)
                    call nextk(type)
                end do
                tried = .true.
                j = 0
                end if
            end if
            omp = mp

            ! calculate next k polynomial and new u and v
            call calcsc(type)
            call nextk(type)
            call calcsc(type)
            call newest(type, ui, vi)

            ! if vi is zero the iteration is not converging
            if (vi == 0.0_wp) exit
            relstp = abs((vi - v)/vi)
            u = ui
            v = vi

        end do main

      end subroutine quadit

      subroutine realit(sss, nz, iflag)

         !! variable-shift h polynomial iteration for a real zero.

         real(wp), intent(inout) :: sss !! starting iterate
         integer, intent(out)    :: nz !! number of zero found
         integer, intent(out)    :: iflag !! flag to indicate a pair of zeros near real axis.

         real(wp) :: pv, kv, t, s, ms, mp, omp, ee
         integer  :: i, j

         nz = 0
         s = sss
         iflag = 0
         j = 0

         ! main loop
         main : do
             pv = p(1)

            ! evaluate p at s
            qp(1) = pv
            do i = 2, nn
                pv = pv*s + p(i)
                qp(i) = pv
            end do
            mp = abs(pv)

            ! compute a rigorous bound on the error in evaluating p
            ms = abs(s)
            ee = (mre/(are + mre))*abs(qp(1))
            do i = 2, nn
                ee = ee*ms + abs(qp(i))
            end do

            ! iteration has converged sufficiently if the
            ! polynomial value is less than 20 times this bound
            if (mp <= 20.0_wp*((are + mre)*ee - mre*mp)) then
                nz = 1
                szr = s
                szi = 0.0_wp
                return
            end if
            j = j + 1

            ! stop iteration after 10 steps
            if (j > 10) return
            if (j >= 2) then
                if (abs(t) <= 0.001_wp * abs(s - t) .and. mp > omp) then
                    ! a cluster of zeros near the real axis has been encountered,
                    ! return with iflag set to initiate a quadratic iteration
                    iflag = 1
                    sss = s
                    return
                end if
            end if

            ! return if the polynomial value has increased significantly
            omp = mp

            ! compute t, the next polynomial, and the new iterate
            kv = k(1)
            qk(1) = kv
            do i = 2, n
                kv = kv*s + k(i)
                qk(i) = kv
            end do
            if (abs(kv) > abs(k(n))*10.0_wp*eta) then
                ! use the scaled form of the recurrence if the value of k at s is nonzero
                t = -pv/kv
                k(1) = qp(1)
                do i = 2, n
                    k(i) = t*qk(i - 1) + qp(i)
                end do
            else
                ! use unscaled form
                k(1) = 0.0_wp
                do i = 2, n
                    k(i) = qk(i - 1)
                end do
            end if
            kv = k(1)
            do i = 2, n
                kv = kv*s + k(i)
            end do
            t = 0.0_wp
            if (abs(kv) > abs(k(n))*10.*eta) t = -pv/kv
            s = s + t

         end do main

      end subroutine realit

      subroutine calcsc(type)

         !! this routine calculates scalar quantities used to
         !! compute the next k polynomial and new estimates of
         !! the quadratic coefficients.

         integer, intent(out)  :: type !! integer variable set here indicating how the
                                       !! calculations are normalized to avoid overflow

         ! synthetic division of k by the quadratic 1,u,v
         call quadsd(n, u, v, k, qk, c, d)
         if (abs(c) <= abs(k(n))*100.0_wp*eta) then
            if (abs(d) <= abs(k(n - 1))*100.0_wp*eta) then
               type = 3
               ! type=3 indicates the quadratic is almost a factor of k
               return
            end if
         end if

         if (abs(d) >= abs(c)) then
            type = 2
            ! type=2 indicates that all formulas are divided by d
            e = a/d
            f = c/d
            g = u*b
            h = v*b
            a3 = (a + g)*e + h*(b/d)
            a1 = b*f - a
            a7 = (f + u)*a + h
            return
         end if
         type = 1
         ! type=1 indicates that all formulas are divided by c
         e = a/c
         f = d/c
         g = u*e
         h = v*b
         a3 = a*e + (h/c + g)*b
         a1 = b - a*(d/c)
         a7 = a + g*d + h*f

      end subroutine calcsc

      subroutine nextk(type)

         !! computes the next k polynomials using scalars computed in calcsc.

         integer, intent(in)  :: type

         real(wp) :: temp
         integer  :: i

         if (type /= 3) then
            temp = a
            if (type == 1) temp = b
            if (abs(a1) <= abs(temp)*eta*10.0_wp) then
               ! if a1 is nearly zero then use a special form of the recurrence
               k(1) = 0.0_wp
               k(2) = -a7*qp(1)
               do i = 3, n
                  k(i) = a3*qk(i - 2) - a7*qp(i - 1)
               end do
               return
            end if

            ! use scaled form of the recurrence
            a7 = a7/a1
            a3 = a3/a1
            k(1) = qp(1)
            k(2) = qp(2) - a7*qp(1)
            do i = 3, n
               k(i) = a3*qk(i - 2) - a7*qp(i - 1) + qp(i)
            end do
            return
         end if

         ! use unscaled form of the recurrence if type is 3
         k(1) = 0.0_wp
         k(2) = 0.0_wp
         do i = 3, n
            k(i) = qk(i - 2)
         end do

      end subroutine nextk

      subroutine newest(type, uu, vv)

         ! compute new estimates of the quadratic coefficients
         ! using the scalars computed in calcsc.

         integer, intent(in)     :: type
         real(wp), intent(out)  :: uu
         real(wp), intent(out)  :: vv

         real(wp) :: a4, a5, b1, b2, c1, c2, c3, c4, temp

         ! use formulas appropriate to setting of type.
         if (type /= 3) then
            if (type /= 2) then
               a4 = a + u*b + h*f
               a5 = c + (u + v*f)*d
            else
               a4 = (a + g)*f + h
               a5 = (f + u)*c + v*d
            end if

            ! evaluate new quadratic coefficients.
            b1 = -k(n)/p(nn)
            b2 = -(k(n - 1) + b1*p(n))/p(nn)
            c1 = v*b2*a1
            c2 = b1*a7
            c3 = b1*b1*a3
            c4 = c1 - c2 - c3
            temp = a5 + b1*a4 - c4
            if (temp /= 0.0_wp) then
               uu = u - (u*(c3 + c2) + v*(b1*a1 + b2*a7))/temp
               vv = v*(1.0_wp+c4/temp)
               return
            end if
         end if

         ! if type=3 the quadratic is zeroed
         uu = 0.0_wp
         vv = 0.0_wp

      end subroutine newest

      subroutine quadsd(nn, u, v, p, q, a, b)

         ! divides `p` by the quadratic `1,u,v` placing the
         ! quotient in `q` and the remainder in `a,b`.

         integer, intent(in)    :: nn
         real(wp), intent(in)   :: u, v, p(nn)
         real(wp), intent(out)  :: q(nn), a, b

         real(wp)  :: c
         integer   :: i

         b = p(1)
         q(1) = b
         a = p(2) - u*b
         q(2) = a
         do i = 3, nn
            c = p(i) - u*a - v*b
            q(i) = c
            b = a
            a = c
         end do

      end subroutine quadsd

      subroutine quad(a, b1, c, sr, si, lr, li)

         !! calculate the zeros of the quadratic a*z**2+b1*z+c.
         !! the quadratic formula, modified to avoid overflow, is used to find the
         !! larger zero if the zeros are real and both zeros are complex.
         !! the smaller real zero is found directly from the product of the zeros c/a.

         real(wp), intent(in)   :: a, b1, c
         real(wp), intent(out)  :: sr, si, lr, li

         real(wp) :: b, d, e

         if (a == 0.0_wp) then
            sr = 0.0_wp
            if (b1 /= 0.0_wp) sr = -c/b1
            lr = 0.0_wp
            si = 0.0_wp
            li = 0.0_wp
            return
         end if

         if (c == 0.0_wp) then
            sr = 0.0_wp
            lr = -b1/a
            si = 0.0_wp
            li = 0.0_wp
            return
         end if

         ! compute discriminant avoiding overflow
         b = b1/2.0_wp
         if (abs(b) >= abs(c)) then
            e = 1.0_wp - (a/b)*(c/b)
            d = sqrt(abs(e))*abs(b)
         else
            e = a
            if (c < 0.0_wp) e = -a
            e = b*(b/abs(c)) - e
            d = sqrt(abs(e))*sqrt(abs(c))
         end if

         if (e >= 0.0_wp) then
            ! real zeros
            if (b >= 0.0_wp) d = -d
            lr = (-b + d)/a
            sr = 0.0_wp
            if (lr /= 0.0_wp) sr = (c/lr)/a
            si = 0.0_wp
            li = 0.0_wp
            return
         end if

         ! complex conjugate zeros
         sr = -b/a
         lr = sr
         si = abs(d/a)
         li = -si

      end subroutine quad

   end subroutine rpoly
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes the roots of a cubic polynomial of the form
!  `a(1) + a(2)*z + a(3)*z**2 + a(4)*z**3 = 0`
!
!### See also
! * A. H. Morris, "NSWC Library of Mathematical Subroutines",
!   Naval Surface Warfare Center, NSWCDD/TR-92/425, January 1993
!
!### History
! * Original version by Alfred H. Morris & William Davis, Naval Surface Weapons Center

   subroutine dcbcrt(a, zr, zi)

      implicit none

      real(wp), dimension(4), intent(in)    :: a    !! coefficients
      real(wp), dimension(3), intent(out)   :: zr   !! real components of roots
      real(wp), dimension(3), intent(out)   :: zi   !! imaginary components of roots

      real(wp) :: arg, c, cf, d, p, p1, q, q1, &
                  r, ra, rb, rq, rt, r1, s, sf, sq, sum, &
                  t, tol, t1, w, w1, w2, x, x1, x2, x3, y, &
                  y1, y2, y3

      real(wp), parameter :: sqrt3 = sqrt(3.0_wp)

      if (a(4) == 0.0_wp) then  ! quadratic equation

         call dqdcrt(a(1:3), zr(1:2), zi(1:2))

         !there are only two roots, so just duplicate the second one:
         zr(3) = zr(2)
         zi(3) = zi(2)

      else

         if (a(1) == 0.0_wp) then ! quadratic

            zr(1) = 0.0_wp
            zi(1) = 0.0_wp
            call dqdcrt(a(2:4), zr(2:3), zi(2:3))

         else

            p = a(3)/(3.0_wp*a(4))
            q = a(2)/a(4)
            r = a(1)/a(4)
            tol = 4.0_wp*eps

            c = 0.0_wp
            t = a(2) - p*a(3)
            if (abs(t) > tol*abs(a(2))) c = t/a(4)

            t = 2.0_wp*p*p - q
            if (abs(t) <= tol*abs(q)) t = 0.0_wp
            d = r + p*t

            if (abs(d) <= tol*abs(r)) then

               !case when d = 0

               zr(1) = -p
               zi(1) = 0.0_wp
               w = sqrt(abs(c))
               if (c < 0.0_wp) then

                  if (p /= 0.0_wp) then

                     x = -(p + sign(w, p))
                     zr(3) = x
                     zi(2) = 0.0_wp
                     zi(3) = 0.0_wp
                     t = 3.0_wp*a(1)/(a(3)*x)

                     if (abs(p) > abs(t)) then
                        zr(2) = zr(1)
                        zr(1) = t
                     else
                        zr(2) = t
                     end if

                  else

                     zr(2) = w
                     zr(3) = -w
                     zi(2) = 0.0_wp
                     zi(3) = 0.0_wp

                  end if

               else

                  zr(2) = -p
                  zr(3) = zr(2)
                  zi(2) = w
                  zi(3) = -w

               end if

            else

               !set  sq = (a(4)/s)**2 * (c**3/27 + d**2/4)

               s = max(abs(a(1)), abs(a(2)), abs(a(3)))
               p1 = a(3)/(3.0_wp*s)
               q1 = a(2)/s
               r1 = a(1)/s

               t1 = q - 2.25_wp*p*p
               if (abs(t1) <= tol*abs(q)) t1 = 0.0_wp
               w = 0.25_wp*r1*r1
               w1 = 0.5_wp*p1*r1*t
               w2 = q1*q1*t1/27.0_wp

               if (w1 < 0.0_wp) then

                  if (w2 < 0.0_wp) then
                     sq = w + (w1 + w2)
                  else
                     w = w + w2
                     sq = w + w1
                  end if

               else

                  w = w + w1
                  sq = w + w2

               end if

               if (abs(sq) <= tol*w) sq = 0.0_wp
               rq = abs(s/a(4))*sqrt(abs(sq))

               if (sq < 0.0_wp) then

                  !all roots are real

                  arg = atan2(rq, -0.5_wp*d)
                  cf = cos(arg/3.0_wp)
                  sf = sin(arg/3.0_wp)
                  rt = sqrt(-c/3.0_wp)
                  y1 = 2.0_wp*rt*cf
                  y2 = -rt*(cf + sqrt3*sf)
                  y3 = -(d/y1)/y2

                  x1 = y1 - p
                  x2 = y2 - p
                  x3 = y3 - p
                  if (abs(x1) > abs(x2)) then
                     t = x1
                     x1 = x2
                     x2 = t
                  end if
                  if (abs(x2) > abs(x3)) then
                     t = x2
                     x2 = x3
                     x3 = t
                     if (abs(x1) > abs(x2)) then
                        t = x1
                        x1 = x2
                        x2 = t
                     end if
                  end if

                  w = x3
                  if (abs(x2) < 0.1_wp*abs(x3)) then
                     call roundoff()
                  else
                     if (abs(x1) < 0.1_wp*abs(x2)) x1 = -(r/x3)/x2
                     zr(1) = x1
                     zr(2) = x2
                     zr(3) = x3
                     zi(1) = 0.0_wp
                     zi(2) = 0.0_wp
                     zi(3) = 0.0_wp
                  end if

               else

                  !real and complex roots

                  ra = dcbrt(-0.5_wp*d - sign(rq, d))
                  rb = -c/(3.0_wp*ra)
                  t = ra + rb
                  w = -p
                  x = -p
                  if (abs(t) > tol*abs(ra)) then
                     w = t - p
                     x = -0.5_wp*t - p
                     if (abs(x) <= tol*abs(p)) x = 0.0_wp
                  end if
                  t = abs(ra - rb)
                  y = 0.5_wp*sqrt3*t

                  if (t > tol*abs(ra)) then

                     if (abs(x) < abs(y)) then
                        s = abs(y)
                        t = x/y
                     else
                        s = abs(x)
                        t = y/x
                     end if

                     if (s < 0.1_wp*abs(w)) then
                        call roundoff()
                     else
                        w1 = w/s
                        sum = 1.0_wp + t*t
                        if (w1*w1 < 1.0e-2_wp*sum) w = -((r/sum)/s)/s
                        zr(1) = w
                        zr(2) = x
                        zr(3) = x
                        zi(1) = 0.0_wp
                        zi(2) = y
                        zi(3) = -y
                     end if

                  else

                     !at least two roots are equal

                     zi(1) = 0.0_wp
                     zi(2) = 0.0_wp
                     zi(3) = 0.0_wp

                     if (abs(x) < abs(w)) then
                        if (abs(x) < 0.1_wp*abs(w)) then
                           call roundoff()
                        else
                           zr(1) = x
                           zr(2) = x
                           zr(3) = w
                        end if
                     else
                        if (abs(w) < 0.1_wp*abs(x)) w = -(r/x)/x
                        zr(1) = w
                        zr(2) = x
                        zr(3) = x
                     end if

                  end if

               end if

            end if

         end if

      end if

   contains

      !*************************************************************
      subroutine roundoff()

        !!  here `w` is much larger in magnitude than the other roots.
        !!  as a result, the other roots may be exceedingly inaccurate
        !!  because of roundoff error. to deal with this, a quadratic
        !!  is formed whose roots are the same as the smaller roots of
        !!  the cubic. this quadratic is then solved.
        !!
        !!  this code was written by william l. davis (nswc).

         implicit none

         real(wp), dimension(3) :: aq

         aq(1) = a(1)
         aq(2) = a(2) + a(1)/w
         aq(3) = -a(4)*w
         call dqdcrt(aq, zr(1:2), zi(1:2))
         zr(3) = w
         zi(3) = 0.0_wp

         if (zi(1) /= 0.0_wp) then
            zr(3) = zr(2)
            zi(3) = zi(2)
            zr(2) = zr(1)
            zi(2) = zi(1)
            zr(1) = w
            zi(1) = 0.0_wp
         end if

      end subroutine roundoff
      !*************************************************************

   end subroutine dcbcrt
!*****************************************************************************************

!*****************************************************************************************
!>
!  Cube root of a real number.

   pure real(wp) function dcbrt(x)

      implicit none

      real(wp), intent(in) :: x

      real(wp), parameter :: third = 1.0_wp/3.0_wp

      dcbrt = sign(abs(x)**third, x)

   end function dcbrt
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes the roots of a quadratic polynomial of the form
!  `a(1) + a(2)*z + a(3)*z**2 = 0`
!
!### See also
! * A. H. Morris, "NSWC Library of Mathematical Subroutines",
!   Naval Surface Warfare Center, NSWCDD/TR-92/425, January 1993

   pure subroutine dqdcrt(a, zr, zi)

      implicit none

      real(wp), dimension(3), intent(in)    :: a    !! coefficients
      real(wp), dimension(2), intent(out)   :: zr   !! real components of roots
      real(wp), dimension(2), intent(out)   :: zi   !! imaginary components of roots

      real(wp) :: d, r, w

      if (a(3) == 0.0_wp) then

         !it is really a linear equation:

         if (a(2) == 0.0_wp) then  !degenerate case, just return zeros

            zr = 0.0_wp
            zi = 0.0_wp

         else

            !there is only one root (real), so just duplicate it:

            zr(1) = -a(1)/a(2)
            zr(2) = zr(1)

            zi = 0.0_wp

         end if

      else

         if (a(1) == 0.0_wp) then

            zr(1) = 0.0_wp
            zi(1) = 0.0_wp
            zr(2) = -a(2)/a(3)
            zi(2) = 0.0_wp

         else

            d = a(2)*a(2) - 4.0_wp*a(1)*a(3)

            if (abs(d) <= 2.0_wp*eps*a(2)*a(2)) then

               !equal real roots

               zr(1) = -0.5_wp*a(2)/a(3)
               zr(2) = zr(1)
               zi(1) = 0.0_wp
               zi(2) = 0.0_wp

            else

               r = sqrt(abs(d))
               if (d < 0.0_wp) then

                  !complex roots

                  zr(1) = -0.5_wp*a(2)/a(3)
                  zr(2) = zr(1)
                  zi(1) = abs(0.5_wp*r/a(3))
                  zi(2) = -zi(1)

               else

                  !distinct real roots

                  zi(1) = 0.0_wp
                  zi(2) = 0.0_wp

                  if (a(2) /= 0.0_wp) then

                     w = -(a(2) + sign(r, a(2)))
                     zr(1) = 2.0_wp*a(1)/w
                     zr(2) = 0.5_wp*w/a(3)

                  else

                     zr(1) = abs(0.5_wp*r/a(3))
                     zr(2) = -zr(1)

                  end if

               end if

            end if

         end if

      end if

   end subroutine dqdcrt
!*****************************************************************************************

!*****************************************************************************************
end module polyroots_module
!*****************************************************************************************