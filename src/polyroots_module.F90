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

   real(wp), parameter :: eps = epsilon(1.0_wp) !! machine epsilon
   real(wp), parameter :: pi = acos(-1.0_wp)

   ! general polynomial routines:
   public :: polyroots
   public :: rpoly
   public :: cpzero
   public :: rpzero
   public :: rpqr79
   public :: cpqr79
   public :: qr_algeq_solver

   ! special polynomial routines:
   public :: dcbcrt
   public :: dqdcrt

   ! utility routines:
   public :: dcbrt
   public :: cpevl
   public :: newton_root_polish

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

    subroutine rpoly(op, degree, zeror, zeroi, istat)

      implicit none

      real(wp), intent(in)    :: op(:) !! vector of coefficients in order of decreasing powers
      integer, intent(in)     :: degree !! degree of polynomial
      real(wp), intent(out)   :: zeror(:) !! output vector of real parts of the zeros
      real(wp), intent(out)   :: zeroi(:) !! output vector of imaginary parts of the zeros
      integer, intent(out)    :: istat !! status output:
                                       !!
                                       !! * 0 : success
                                       !! * -1 : leading coefficient is zero
                                       !! * -2 : no roots found
                                       !! * >0 : the number of zeros found

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

      real(wp),parameter :: deg2rad = pi / 180.0_wp
      real(wp),parameter :: cosr = cos(94.0_wp * deg2rad)
      real(wp),parameter :: sinr = sin(86.0_wp * deg2rad)
      real(wp),parameter :: base = radix(0.0_wp)
      real(wp),parameter :: eta = eps
      real(wp),parameter :: infin = huge(0.0_wp)
      real(wp),parameter :: smalno = tiny(0.0_wp)
      real(wp),parameter :: sqrthalf = sqrt(0.5_wp)
      real(wp),parameter :: are = eta !! unit error in +
      real(wp),parameter :: mre = eta !! unit error in *
      real(wp),parameter :: lo = smalno/eta

      ! initialization of constants for shift rotation
      xx = sqrthalf
      yy = -xx
      istat = 0
      n = degree
      nn = n + 1

      ! algorithm fails if the leading coefficient is zero.
      if (op(1) == 0.0_wp) then
         istat = -1
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
      istat = degree - n
      if (istat==0) istat = -2  ! if not roots found

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
         logical  :: vpass, spass, vtry, stry, skip

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
                  skip = (spass .and. ((.not. vpass) .or. tss < tvv))

                  do
                    do

                      if (.not. skip) then
                        call quadit(ui, vi, nz)
                        if (nz > 0) return

                        ! quadratic iteration has failed. flag that it has
                        ! been tried and decrease the convergence criterion.
                        vtry = .true.
                        betav = betav * 0.25_wp

                        ! try linear iteration if it has not been tried and
                        ! the s sequence is converging
                        if (stry .or. (.not. spass)) exit
                        k(1:n) = svk(1:n)
                      end if

                      call realit(s, nz, iflag)
                      if (nz > 0) return

                      ! linear iteration has failed.  flag that it has been
                      ! tried and decrease the convergence criterion
                      stry = .true.
                      betas = betas * 0.25_wp
                      if (iflag == 0) exit

                      ! if linear iteration signals an almost double real
                      ! zero attempt quadratic interation
                      ui = -(s + s)
                      vi = s*s

                    end do

                    ! restore variables
                    u = svu
                    v = svv
                    k(1:n) = svk(1:n)

                    ! try quadratic iteration if it has not been tried
                    ! and the v sequence is converging
                    if (.not.(vpass .and. (.not. vtry))) exit

                  end do

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
!>
!  Solve the real coefficient algebraic equation by the qr-method.
!
!### Reference
!  * [`/opt/companion.tgz`](https://netlib.org/opt/companion.tgz) on Netlib
!    from [Edelman & Murakami (1995)](https://www.ams.org/journals/mcom/1995-64-210/S0025-5718-1995-1262279-2/S0025-5718-1995-1262279-2.pdf),

   subroutine qr_algeq_solver(n,c,zr,zi,detil,istatus)

    implicit none

    integer,intent(in) :: n !! degree of the monic polynomial.
    real(wp),intent(in) :: c(n+1) !! coefficients of the polynomial. in order of decreasing powers.
    real(wp),intent(out) :: zr(n) !! real part of output roots
    real(wp),intent(out) :: zi(n) !! imaginary part of output roots
    real(wp),intent(out) :: detil !! accuracy hint.
    integer,intent(out) :: istatus !! return code:
                                   !!
                                   !! * -1 : degree <= 0
                                   !! * -2 : leading coefficient `c(1)` is zero
                                   !! * 0 : success
                                   !! * otherwise, the return code from `hqr_eigen_hessenberg`

    real(wp),allocatable :: a(:,:) !! work matrix
    integer ,allocatable :: cnt(:) !! work area for counting the qr-iterations
    integer :: i
    integer :: iter
    real(wp) :: afnorm

    ! check for invalid arguments
    if ( n<=0 ) then
       istatus = -1
       return
    endif
    if (c(1) == 0.0_wp) then
        ! leading coefficient is zero.
        istatus = -2
        return
    end if

    allocate(a(n,n))
    allocate(cnt(n))

    ! build the companion matrix a.
    call build_companion(n,a,c)

    ! balancing the a itself.
    call balance_companion(n,a)

    ! compute the frobenius norm of the balanced companion matrix a.
    afnorm = frobenius_norm_companion(n,a)

    ! qr iterations from a.
    call hqr_eigen_hessenberg(n,a,zr,zi,cnt,istatus)
    if ( istatus/=0 ) then
        write(*,'(A,1X,I4)')  'abnormal from hqr_eigen_hessenberg. istatus=' , istatus
        if ( istatus==1 ) write(*,'(A)') 'matrix is completely zero.'
        if ( istatus==2 ) write(*,'(A)') 'qr iteration does not converged.'
        if ( istatus>3 )  write(*,'(A)') 'arguments violate the condition.'
        return
    endif

    ! count the total qr iteration.
    iter = 0
    do i = 1 , n
       if ( cnt(i)>0 ) iter = iter + cnt(i)
    enddo

    ! calculate the accuracy hint.
    detil = eps*n*iter*afnorm

contains

subroutine build_companion(n,a,c)
    !!  build the companion matrix of the polynomial.
    !!  (this was modified to allow for non-monic polynomials)
    implicit none

    integer,intent(in)   :: n
    real(wp),intent(out) :: a(n,n)
    real(wp),intent(in)  :: c(n+1) !! coefficients in order of decreasing powers

    integer :: i !! counter

    ! create the companion matrix
    a = 0.0_wp
    do i = 1, n-1
        a(i+1,i) = 1.0_wp
    end do
    do i = n,1,-1
        a(n-i+1,n) = -c(i+1) / c(1)
    end do

end subroutine build_companion

subroutine balance_companion(n,a)

    !!  blancing the unsymmetric matrix `a`.
    !!
    !!  this fortran code is based on the algol code "balance" from paper:
    !!   "balancing a matrixfor calculation of eigenvalues and eigenvectors"
    !!   by b.n.parlett and c.reinsch, numer. math. 13, 293-304(1969).
    !!
    !!  note: the only non-zero elements of the companion matrix are touched.

    implicit none

    integer,intent(in) :: n
    real(wp),intent(inout) :: a(n,n)

    integer,parameter :: b = radix(1.0_wp) !! base of the floating point representation on the machine
    integer,parameter :: b2 = b**2

    integer :: i , j
    real(wp) :: c , f , g , r , s
    logical :: noconv

    if ( n<=1 ) return ! do nothing

    ! iteration:
    main : do
        noconv = .false.
        do i = 1 , n
            ! touch only non-zero elements of companion.
            if ( i/=n ) then
                c = abs(a(i+1,i))
            else
                c = 0.0_wp
                do j = 1 , n - 1
                    c = c + abs(a(j,n))
                enddo
            endif
            if ( i==1 ) then
                r = abs(a(1,n))
            elseif ( i/=n ) then
                r = abs(a(i,i-1)) + abs(a(i,n))
            else
                r = abs(a(i,i-1))
            endif

            if ( c/=0.0_wp .and. r/=0.0_wp ) then

                g = r/b
                f = 1.0_wp
                s = c + r
                do
                    if ( c>=g ) exit
                    f = f*b
                    c = c*b2
                end do
                g = r*b
                do
                    if ( c<g ) exit
                    f = f/b
                    c = c/b2
                end do
                if ( (c+r)/f<0.95_wp*s ) then
                    g = 1.0_wp/f
                    noconv = .true.
                    ! touch only non-zero elements of companion.
                    if ( i==1 ) then
                        a(1,n) = a(1,n)*g
                    else
                        a(i,i-1) = a(i,i-1)*g
                        a(i,n) = a(i,n)*g
                    endif
                    if ( i/=n ) then
                        a(i+1,i) = a(i+1,i)*f
                    else
                        do j = 1 , n
                            a(j,i) = a(j,i)*f
                        enddo
                    endif
                endif
            endif
        enddo
        if ( noconv ) cycle main
        exit main
    end do main

    end subroutine balance_companion

    function frobenius_norm_companion(n,a) result(afnorm)

    !!  calculate the frobenius norm of the companion-like matrix.
    !!  note: the only non-zero elements of the companion matrix are touched.

    implicit none

    integer,intent(in) :: n
    real(wp),intent(in) :: a(n,n)
    real(wp) :: afnorm

    integer :: i , j
    real(wp) :: sum

    sum = 0.0_wp
    do j = 1 , n - 1
       sum = sum + a(j+1,j)**2
    enddo
    do i = 1 , n
       sum = sum + a(i,n)**2
    enddo
    afnorm = sqrt(sum)

    end function frobenius_norm_companion

    subroutine hqr_eigen_hessenberg(n0,h,wr,wi,cnt,istatus)

    !!  eigenvalue computation by the householder qr method
    !!  for the real hessenberg matrix.
    !!
    !! this fortran code is based on the algol code "hqr" from the paper:
    !!       "the qr algorithm for real hessenberg matrices"
    !!       by r.s.martin, g.peters and j.h.wilkinson,
    !!       numer. math. 14, 219-231(1970).
    !!
    !! comment: finds the eigenvalues of a real upper hessenberg matrix,
    !!          h, stored in the array h(1:n0,1:n0), and stores the real
    !!          parts in the array wr(1:n0) and the imaginary parts in the
    !!          array wi(1:n0).
    !!          the procedure fails if any eigenvalue takes more than
    !!          30 iterations.

    implicit none

    integer,intent(in) :: n0
    real(wp) :: h(n0,n0)
    real(wp) :: wr(n0)
    real(wp) :: wi(n0)
    integer :: cnt(n0)
    integer,intent(out) :: istatus

    integer :: i , j , k , l , m , na , its
    real(wp) :: p , q , r , s , t , w , x , y , z
    logical :: notlast
    integer :: n

    ! note: n is changing in this subroutine.
    n = n0
    istatus = 0
    t = 0.0_wp

100 if ( n==0 ) return

    its = 0
    na = n - 1
    ! look for single small sub-diagonal element
200 do l = n , 2 , -1
       if ( abs(h(l,l-1))<=eps*(abs(h(l-1,l-1))+abs(h(l,l))) ) goto 300
    enddo
    l = 1

300 x = h(n,n)
    if ( l==n ) then
       ! one root found
       wr(n) = x + t
       wi(n) = 0.0_wp
       cnt(n) = its
       n = na
       goto 100
    else
       y = h(na,na)
       w = h(n,na)*h(na,n)
       if ( l==na ) then
          ! comment: two roots found
          p = (y-x)/2
          q = p**2 + w
          y = sqrt(abs(q))
          cnt(n) = -its
          cnt(na) = its
          x = x + t
          if ( q>0.0_wp ) then
             ! real pair
             if ( p<0.0_wp ) y = -y
             y = p + y
             wr(na) = x + y
             wr(n) = x - w/y
             wi(na) = 0.0_wp
             wi(n) = 0.0_wp
          else
            ! complex pair
             wr(na) = x + p
             wr(n) = x + p
             wi(na) = y
             wi(n) = -y
          endif
          n = n - 2
          goto 100
       else
          if ( its==30 ) then
             istatus = 1
             return
          endif
          if ( its==10 .or. its==20 ) then
             ! form exceptional shift
             t = t + x
             do i = 1 , n
                h(i,i) = h(i,i) - x
             enddo
             s = abs(h(n,na)) + abs(h(na,n-2))
             y = 0.75_wp*s
             x = y
             w = -0.4375_wp*s**2
          endif
          its = its + 1
          ! look for two consecutive small sub-diagonal elements
          do m = n - 2 , l , -1
             z = h(m,m)
             r = x - z
             s = y - z
             p = (r*s-w)/h(m+1,m) + h(m,m+1)
             q = h(m+1,m+1) - z - r - s
             r = h(m+2,m+1)
             s = abs(p) + abs(q) + abs(r)
             p = p/s
             q = q/s
             r = r/s
             if ( m==l ) exit
             if ( abs(h(m,m-1))*(abs(q)+abs(r))<=eps*abs(p) &
                  *(abs(h(m-1,m-1))+abs(z)+abs(h(m+1,m+1))) ) exit
          enddo

          do i = m + 2 , n
             h(i,i-2) = 0.0_wp
          enddo
          do i = m + 3 , n
             h(i,i-3) = 0.0_wp
          enddo
          ! double qr step involving rows l to n and columns m to n
          do k = m , na
             notlast = (k/=na)
             if ( k/=m ) then
                p = h(k,k-1)
                q = h(k+1,k-1)
                if ( notlast ) then
                   r = h(k+2,k-1)
                else
                   r = 0.0_wp
                endif
                x = abs(p) + abs(q) + abs(r)
                if ( x==0.0_wp ) cycle
                p = p/x
                q = q/x
                r = r/x
             endif
             s = sqrt(p**2+q**2+r**2)
             if ( p<0.0_wp ) s = -s
             if ( k/=m ) then
                h(k,k-1) = -s*x
             elseif ( l/=m ) then
                h(k,k-1) = -h(k,k-1)
             endif
             p = p + s
             x = p/s
             y = q/s
             z = r/s
             q = q/p
             r = r/p
             ! row modification
             do j = k , n
                p = h(k,j) + q*h(k+1,j)
                if ( notlast ) then
                   p = p + r*h(k+2,j)
                   h(k+2,j) = h(k+2,j) - p*z
                endif
                h(k+1,j) = h(k+1,j) - p*y
                h(k,j) = h(k,j) - p*x
             enddo
             if ( k+3<n ) then
                j = k + 3
             else
                j = n
             endif
             ! column modification;
             do i = l , j
                p = x*h(i,k) + y*h(i,k+1)
                if ( notlast ) then
                   p = p + z*h(i,k+2)
                   h(i,k+2) = h(i,k+2) - p*r
                endif
                h(i,k+1) = h(i,k+1) - p*q
                h(i,k) = h(i,k) - p
             enddo
           enddo
          goto 200
       endif
    endif

    end subroutine hqr_eigen_hessenberg

    end subroutine qr_algeq_solver
!*****************************************************************************************

!*****************************************************************************************
!>
!  Evaluate a complex polynomial and its derivatives.
!  Optionally compute error bounds for these values.
!
!### REVISION HISTORY  (YYMMDD)
!  * 810223  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890831  Modified array declarations.  (WRB)
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900402  Added TYPE section.  (WRB)
!  * Jacob Williams, 9/13/2022 : modernized this routine

    subroutine cpevl(n,m,a,z,c,b,kbd)

    implicit none

    integer,intent(in) :: n !! Degree of the polynomial
    integer,intent(in) :: m !! Number of derivatives to be calculated:
                            !!
                            !!  * M=0 evaluates only the function
                            !!  * M=1 evaluates the function and first derivative, etc.
                            !!
                            !! if M > N+1 function and all N derivatives will be calculated.
    complex(wp),intent(in) :: a(*) !! vector containing the N+1 coefficients of polynomial.
                                   !! A(I) = coefficient of Z**(N+1-I)
    complex(wp),intent(in) :: z !! point at which the evaluation is to take place
    complex(wp),intent(out) :: c(*) !! Array of 2(M+1) words: C(I+1) contains the complex value of the I-th
                                    !! derivative at Z, I=0,...,M
    complex(wp),intent(out) :: b(*) !! Array of 2(M+1) words: B(I) contains the bounds on the real and imaginary parts
                                    !! of C(I) if they were requested. only needed if bounds are to be calculated.
                                    !! It is not used otherwise.
    logical,intent(in) :: kbd !! A logical variable, e.g. .TRUE. or .FALSE. which is
                              !! to be set .TRUE. if bounds are to be computed.

    real(wp) :: r , s
    integer :: i , j , mini , np1
    complex(wp) :: ci , cim1 , bi , bim1 , t , za , q

    real(wp),parameter :: d1 = real(radix(1.0_wp))**(1-digits(1.0_wp))

    za(q) = cmplx(abs(real(q,wp)),abs(aimag(q)),wp)
    np1 = n + 1
    do j = 1 , np1
        ci = 0.0_wp
        cim1 = a(j)
        bi = 0.0_wp
        bim1 = 0.0_wp
        mini = min(m+1,n+2-j)
        do i = 1 , mini
            if ( j/=1 ) ci = c(i)
            if ( i/=1 ) cim1 = c(i-1)
            c(i) = cim1 + z*ci
            if ( kbd ) then
                if ( j/=1 ) bi = b(i)
                if ( i/=1 ) bim1 = b(i-1)
                t = bi + (3.0_wp*d1+4.0_wp*d1*d1)*za(ci)
                r = real(za(z)*cmplx(real(t,wp),-aimag(t), wp), wp)
                s = aimag(za(z)*t)
                b(i) = (1.0_wp + 8.0_wp*d1)*(bim1+d1*za(cim1)+cmplx(r,s,wp))
                if ( j==1 ) b(i) = 0.0_wp
            endif
        enddo
    enddo

    end subroutine cpevl
!*****************************************************************************************

!*****************************************************************************************
!>
!  Find the zeros of a polynomial with complex coefficients:
!  `P(Z)= A(1)*Z**N + A(2)*Z**(N-1) +...+ A(N+1)`
!
!### REVISION HISTORY  (YYMMDD)
!  * 810223  DATE WRITTEN. Kahaner, D. K., (NBS)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890531  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * Jacob Williams, 9/13/2022 : modernized this routine

subroutine cpzero(in,a,r,t,iflg,s)

    implicit none

    integer,intent(in) :: in !! degree of P(Z)
    complex(wp),intent(in) :: a(*) !! complex vector containing coefficients of P(Z),
                                   !! A(I) = coefficient of Z**(N+1-i)
    complex(wp),intent(inout) :: r(*) !! N word complex vector. On input: containing initial estimates for zeros
                                      !! if these are known. On output: Ith zero
    complex(wp) :: t(*) !! 4(N+1) word array used for temporary storage
    integer,intent(inout) :: iflg !!### On Input:
                                  !!
                                  !! flag to indicate if initial estimates of zeros are input:
                                  !!
                                  !!  * If IFLG == 0, no estimates are input.
                                  !!  * If IFLG /= 0, the vector R contains estimates of the zeros
                                  !!
                                  !! ** WARNING ****** If estimates are input, they must
                                  !!                   be separated, that is, distinct or
                                  !!                   not repeated.
                                  !!### On Output:
                                  !!
                                  !! Error Diagnostics:
                                  !!
                                  !! * If IFLG == 0 on return, all is well
                                  !! * If IFLG == 1 on return, A(1)=0.0 or N=0 on input
                                  !! * If IFLG == 2 on return, the program failed to converge
                                  !!   after 25*N iterations.  Best current estimates of the
                                  !!   zeros are in R(I).  Error bounds are not calculated.
    real(wp),intent(out) :: s(*) !! an N word array. S(I) = bound for R(I)

    integer :: i , imax , j , n , n1 , nit , nmax , nr
    real(wp) :: u , v , x
    complex(wp) :: pn , temp
    complex(wp) :: ctmp(1), btmp(1)

    if ( in<=0 .or. abs(a(1))==0.0_wp ) then
        iflg = 1
    else
        ! check for easily obtained zeros
        n = in
        n1 = n + 1
        if ( iflg==0 ) then
20         n1 = n + 1
            if ( n<=1 ) then
                r(1) = -a(2)/a(1)
                s(1) = 0.0_wp
                return
            elseif ( abs(a(n1))/=0.0_wp ) then
                ! if initial estimates for zeros not given, find some
                temp = -a(2)/(a(1)*n)
                call cpevl(n,n,a,temp,t,t,.false.)
                imax = n + 2
                t(n1) = abs(t(n1))
                do i = 2 , n1
                    t(n+i) = -abs(t(n+2-i))
                    if ( real(t(n+i),wp)<real(t(imax),wp) ) imax = n + i
                enddo
                x = (-real(t(imax),wp)/real(t(n1),wp))**(1.0_wp/(imax-n1))
                do
                    x = 2.0_wp*x
                    call cpevl(n,0,t(n1),cmplx(x,0.0_wp,wp),ctmp,btmp,.false.)
                    pn = ctmp(1)
                    if ( real(pn,wp)>=0.0_wp ) exit
                end do
                u = 0.5_wp*x
                v = x
                do
                    x = 0.5_wp*(u+v)
                    call cpevl(n,0,t(n1),cmplx(x,0.0_wp,wp),ctmp,btmp,.false.)
                    pn = ctmp(1)
                    if ( real(pn,wp)>0.0_wp ) v = x
                    if ( real(pn,wp)<=0.0_wp ) u = x
                    if ( (v-u)<=0.001_wp*(1.0_wp+v) ) exit
                end do
                do i = 1 , n
                    u = (pi/n)*(2*i-1.5_wp)
                    r(i) = max(x,0.001_wp*abs(temp))*cmplx(cos(u),sin(u),wp) + temp
                enddo
            else
                r(n) = 0.0_wp
                s(n) = 0.0_wp
                n = n - 1
                goto 20
            endif
        endif

        ! main iteration loop starts here
        nr = 0
        nmax = 25*n
        do nit = 1 , nmax
            do i = 1 , n
                if ( nit==1 .or. abs(t(i))/=0.0_wp ) then
                    call cpevl(n,0,a,r(i),ctmp,btmp,.true.)
                    pn = ctmp(1)
                    temp = btmp(1)
                    if ( abs(real(pn,wp))+abs(aimag(pn))>real(temp,wp)+aimag(temp) ) then
                        temp = a(1)
                        do j = 1 , n
                            if ( j/=i ) temp = temp*(r(i)-r(j))
                        enddo
                        t(i) = pn/temp
                    else
                        t(i) = 0.0_wp
                        nr = nr + 1
                    endif
                endif
            enddo
            do i = 1 , n
                r(i) = r(i) - t(i)
            enddo
            if ( nr==n ) then
                ! calculate error bounds for zeros
                do nr = 1 , n
                    call cpevl(n,n,a,r(nr),t,t(n+2),.true.)
                    x = abs(cmplx(abs(real(t(1),wp)),abs(aimag(t(1))),wp)+t(n+2))
                    s(nr) = 0.0_wp
                    do i = 1 , n
                        x = x*real(n1-i,wp)/i
                        temp = cmplx(max(abs(real(t(i+1),wp))-real(t(n1+i),wp),0.0_wp), &
                                max(abs(aimag(t(i+1)))-aimag(t(n1+i)),0.0_wp), wp)
                        s(nr) = max(s(nr),(abs(temp)/x)**(1.0_wp/i))
                    enddo
                    s(nr) = 1.0_wp/s(nr)
                enddo
                return
            end if
        enddo
        iflg = 2  ! error exit
    endif

    end subroutine cpzero
!*****************************************************************************************

!*****************************************************************************************
!>
!  Find the zeros of a polynomial with real coefficients:
!  `P(X)= A(1)*X**N + A(2)*X**(N-1) +...+ A(N+1)`
!
!### REVISION HISTORY  (YYMMDD)
!  * 810223  DATE WRITTEN. Kahaner, D. K., (NBS)
!  * 890206  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * Jacob Williams, 9/13/2022 : modernized this routine
!
!@note This is just a wrapper to [[cpzero]]

    subroutine rpzero(n,a,r,t,iflg,s)

    implicit none

    integer,intent(in) :: n !! degree of P(X)
    real(wp),intent(in) :: a(*) !! real vector containing coefficients of P(X),
                                !! A(I) = coefficient of X**(N+1-I)
    complex(wp),intent(inout) :: r(*) !! N word complex vector. On Input: containing initial estimates for zeros
                                      !! if these are known. On output: ith zero.
    complex(wp) :: t(*) !! 6(N+1) word array used for temporary storage
    integer,intent(inout) :: iflg !!### On Input:
                                  !!
                                  !! flag to indicate if initial estimates of zeros are input:
                                  !!
                                  !!  * If IFLG == 0, no estimates are input.
                                  !!  * If IFLG /= 0, the vector R contains estimates of the zeros
                                  !!
                                  !! ** WARNING ****** If estimates are input, they must
                                  !!                   be separated, that is, distinct or
                                  !!                   not repeated.
                                  !!### On Output:
                                  !!
                                  !! Error Diagnostics:
                                  !!
                                  !! * If IFLG == 0 on return, all is well
                                  !! * If IFLG == 1 on return, A(1)=0.0 or N=0 on input
                                  !! * If IFLG == 2 on return, the program failed to converge
                                  !!   after 25*N iterations.  Best current estimates of the
                                  !!   zeros are in R(I).  Error bounds are not calculated.
    real(wp),intent(out) :: s(*) !! an N word array. bound for R(I).

    integer :: i, n1

    n1 = n + 1
    do i = 1 , n1
        t(i) = cmplx(a(i), 0.0_wp, wp)
    enddo
    call cpzero(n,t,r,t(n+2),iflg,s)

    end subroutine rpzero
!*****************************************************************************************

!*****************************************************************************************
!>
!  This routine computes all zeros of a polynomial of degree NDEG
!  with real coefficients by computing the eigenvalues of the
!  companion matrix.
!
!### REVISION HISTORY  (YYMMDD)
!
!  * 800601  DATE WRITTEN. Vandevender, W. H., (SNLA)
!  * 890505  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 911010  Code reworked and simplified.  (RWC and WRB)
!  * Jacob Williams, 9/13/2022 : modernized this routine

    subroutine rpqr79(ndeg,coeff,root,ierr)

    implicit none

    integer,intent(in) :: ndeg !! degree of polynomial
    real(wp),intent(in) :: coeff(*) !! `NDEG+1` coefficients in descending order.  i.e.,
                                    !! `P(Z) = COEFF(1)*(Z**NDEG) + COEFF(NDEG)*Z + COEFF(NDEG+1)`
    complex(wp),intent(out) :: root(*) !! `NDEG` vector of roots
    integer,intent(out) :: ierr !! Output Error Code
                                !!
                                !!### Normal Code:
                                !!
                                !!  * 0 -- means the roots were computed.
                                !!
                                !!### Abnormal Codes
                                !!
                                !!  * 1 -- more than 30 QR iterations on some eigenvalue of the
                                !!    companion matrix
                                !!  * 2 -- COEFF(1)=0.0
                                !!  * 3 -- NDEG is invalid (less than or equal to 0)

    real(wp) :: scale
    integer :: k , kh , kwr , kwi , kcol, km1 , kwend
    real(wp),dimension(:),allocatable :: work !! work array of dimension `NDEG*(NDEG+2)`

    ierr = 0
    if ( abs(coeff(1))==0.0_wp ) then
        ierr = 2
        write(*,*) 'leading coefficient is zero.'
        return
    endif

    if ( ndeg<=0 ) then
        ierr = 3
        write(*,*) 'degree invalid.'
        return
    endif

    if ( ndeg==1 ) then
        root(1) = cmplx(-coeff(2)/coeff(1),0.0_wp, wp)
        return
    endif

    allocate(work(ndeg*(ndeg+2))) ! work array

    scale = 1.0_wp/coeff(1)
    kh = 1
    kwr = kh + ndeg*ndeg
    kwi = kwr + ndeg
    kwend = kwi + ndeg - 1

    do k = 1 , kwend
        work(k) = 0.0_wp
    enddo

    do k = 1 , ndeg
        kcol = (k-1)*ndeg + 1
        work(kcol) = -coeff(k+1)*scale
        if ( k/=ndeg ) work(kcol+k) = 1.0_wp
    enddo

    call hqr(ndeg,ndeg,1,ndeg,work(kh),work(kwr),work(kwi),ierr)

    if ( ierr/=0 ) then
        ierr = 1
        write(*,*) 'no convergence in 30 qr iterations.'
        return
    endif

    do k = 1 , ndeg
        km1 = k - 1
        root(k) = cmplx(work(kwr+km1),work(kwi+km1),wp)
    enddo

    end subroutine rpqr79
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subroutine finds the eigenvalues of a real
!  upper hessenberg matrix by the qr method.
!
!### Reference
!  * this subroutine is a translation of the algol procedure hqr,
!    num. math. 14, 219-231(1970) by martin, peters, and wilkinson.
!    handbook for auto. comp., vol.ii-linear algebra, 359-371(1971).
!
!### History
!  * this version dated september 1989.
!    RESTORED CORRECT INDICES OF LOOPS (200,210,230,240). (9/29/89 BSG).
!    questions and comments should be directed to burton s. garbow,
!    mathematics and computer science div, argonne national laboratory
!  * Jacob Williams, 9/13/2022 : modernized this routine
!
!@note This routine is from [EISPACK](https://netlib.org/eispack/hqr.f)

    subroutine hqr(nm,n,low,igh,h,wr,wi,ierr)

    implicit none

    integer,intent(in) :: nm !! must be set to the row dimension of two-dimensional
                             !! array parameters as declared in the calling program
                             !! dimension statement.
    integer,intent(in) :: n !! order of the matrix
    integer,intent(in) :: low   !! low and igh are integers determined by the balancing
                                !! subroutine  balanc.  if  balanc  has not been used,
                                !! set low=1, igh=n.
    integer,intent(in) :: igh   !! low and igh are integers determined by the balancing
                                !! subroutine  balanc.  if  balanc  has not been used,
                                !! set low=1, igh=n.
    real(wp),intent(inout) :: h(nm,n) !! On input: contains the upper hessenberg matrix.  information about
                                      !! the transformations used in the reduction to hessenberg
                                      !! form by  elmhes  or  orthes, if performed, is stored
                                      !! in the remaining triangle under the hessenberg matrix.
                                      !!
                                      !! On output: has been destroyed.  therefore, it must be saved
                                      !! before calling `hqr` if subsequent calculation and
                                      !! back transformation of eigenvectors is to be performed.
    real(wp),intent(out) :: wr(n) !! the real parts of the eigenvalues.  the eigenvalues
                                  !! are unordered except that complex conjugate pairs
                                  !! of values appear consecutively with the eigenvalue
                                  !! having the positive imaginary part first.  if an
                                  !! error exit is made, the eigenvalues should be correct
                                  !! for indices ierr+1,...,n.
    real(wp),intent(out) :: wi(n) !! the imaginary parts of the eigenvalues.  the eigenvalues
                                  !! are unordered except that complex conjugate pairs
                                  !! of values appear consecutively with the eigenvalue
                                  !! having the positive imaginary part first.  if an
                                  !! error exit is made, the eigenvalues should be correct
                                  !! for indices ierr+1,...,n.
    integer,intent(out) :: ierr !! is set to:
                                !!
                                !!  * zero -- for normal return,
                                !!  * j -- if the limit of 30*n iterations is exhausted
                                !!    while the j-th eigenvalue is being sought.

    integer :: i , j , k , l , m , en , ll , mm , na , &
               itn , its , mp2 , enm2
    real(wp) :: p , q , r , s , t , w , x , y , zz , norm , &
                tst1 , tst2
    logical :: notlas

    ierr = 0
    norm = 0.0_wp
    k = 1

    ! store roots isolated by balance and compute matrix norm
    do i = 1 , n
        do j = k , n
            norm = norm + abs(h(i,j))
        enddo
        k = i
        if ( i<low .or. i>igh ) then
            wr(i) = h(i,i)
            wi(i) = 0.0_wp
        endif
    enddo

    en = igh
    t = 0.0_wp
    itn = 30*n

    do
        ! search for next eigenvalues
        if ( en<low ) return
        its = 0
        na = en - 1
        enm2 = na - 1
        do
            ! look for single small sub-diagonal element
            ! for l=en step -1 until low do --
            do ll = low , en
                l = en + low - ll
                if ( l==low ) exit
                s = abs(h(l-1,l-1)) + abs(h(l,l))
                if ( s==0.0_wp ) s = norm
                tst1 = s
                tst2 = tst1 + abs(h(l,l-1))
                if ( tst2==tst1 ) exit
            enddo
            ! form shift
            x = h(en,en)
            if ( l==en ) then
                ! one root found
                wr(en) = x + t
                wi(en) = 0.0_wp
                en = na
            else
                y = h(na,na)
                w = h(en,na)*h(na,en)
                if ( l==na ) then
                    ! two roots found
                    p = (y-x)/2.0_wp
                    q = p*p + w
                    zz = sqrt(abs(q))
                    x = x + t
                    if ( q<0.0_wp ) then
                        ! complex pair
                        wr(na) = x + p
                        wr(en) = x + p
                        wi(na) = zz
                        wi(en) = -zz
                    else
                        ! real pair
                        zz = p + sign(zz,p)
                        wr(na) = x + zz
                        wr(en) = wr(na)
                        if ( zz/=0.0_wp ) wr(en) = x - w/zz
                        wi(na) = 0.0_wp
                        wi(en) = 0.0_wp
                    endif
                    en = enm2
                elseif ( itn==0 ) then
                    ! set error -- all eigenvalues have not
                    ! converged after 30*n iterations
                    ierr = en
                    return
                else
                    if ( its==10 .or. its==20 ) then
                        ! form exceptional shift
                        t = t + x
                        do i = low , en
                            h(i,i) = h(i,i) - x
                        enddo
                        s = abs(h(en,na)) + abs(h(na,enm2))
                        x = 0.75_wp*s
                        y = x
                        w = -0.4375_wp*s*s
                    endif
                    its = its + 1
                    itn = itn - 1
                    ! look for two consecutive small
                    ! sub-diagonal elements.
                    ! for m=en-2 step -1 until l do --
                    do mm = l , enm2
                        m = enm2 + l - mm
                        zz = h(m,m)
                        r = x - zz
                        s = y - zz
                        p = (r*s-w)/h(m+1,m) + h(m,m+1)
                        q = h(m+1,m+1) - zz - r - s
                        r = h(m+2,m+1)
                        s = abs(p) + abs(q) + abs(r)
                        p = p/s
                        q = q/s
                        r = r/s
                        if ( m==l ) exit
                        tst1 = abs(p)*(abs(h(m-1,m-1))+abs(zz)+abs(h(m+1,m+1)))
                        tst2 = tst1 + abs(h(m,m-1))*(abs(q)+abs(r))
                        if ( tst2==tst1 ) exit
                    enddo

                    mp2 = m + 2

                    do i = mp2 , en
                        h(i,i-2) = 0.0_wp
                        if ( i/=mp2 ) h(i,i-3) = 0.0_wp
                    enddo
                    ! double qr step involving rows l to en and
                    ! columns m to en
                    do k = m , na
                        notlas = k/=na
                        if ( k/=m ) then
                            p = h(k,k-1)
                            q = h(k+1,k-1)
                            r = 0.0_wp
                            if ( notlas ) r = h(k+2,k-1)
                            x = abs(p) + abs(q) + abs(r)
                            if ( x==0.0_wp ) cycle
                            p = p/x
                            q = q/x
                            r = r/x
                        endif
                        s = sign(sqrt(p*p+q*q+r*r),p)
                        if ( k==m ) then
                            if ( l/=m ) h(k,k-1) = -h(k,k-1)
                        else
                            h(k,k-1) = -s*x
                        endif
                        p = p + s
                        x = p/s
                        y = q/s
                        zz = r/s
                        q = q/p
                        r = r/p
                        if ( notlas ) then
                            ! row modification
                            do j = k , en
                                p = h(k,j) + q*h(k+1,j) + r*h(k+2,j)
                                h(k,j) = h(k,j) - p*x
                                h(k+1,j) = h(k+1,j) - p*y
                                h(k+2,j) = h(k+2,j) - p*zz
                            enddo
                            j = min(en,k+3)
                            ! column modification
                            do i = l , j
                                p = x*h(i,k) + y*h(i,k+1) + zz*h(i,k+2)
                                h(i,k) = h(i,k) - p
                                h(i,k+1) = h(i,k+1) - p*q
                                h(i,k+2) = h(i,k+2) - p*r
                            enddo
                        else
                            ! row modification
                            do j = k , en
                                p = h(k,j) + q*h(k+1,j)
                                h(k,j) = h(k,j) - p*x
                                h(k+1,j) = h(k+1,j) - p*y
                            enddo
                            j = min(en,k+3)
                            ! column modification
                            do i = l , j
                                p = x*h(i,k) + y*h(i,k+1)
                                h(i,k) = h(i,k) - p
                                h(i,k+1) = h(i,k+1) - p*q
                            enddo
                        endif
                    enddo
                    cycle
                endif
            endif
            exit
        end do

    end do

    end subroutine hqr
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve for the roots of a real polynomial equation by
!  computing the eigenvalues of the companion matrix.
!
!  This one uses LAPACK for the eigen solver (`sgeev` or `dgeev`).
!
!### Reference
!  * Code from ivanpribec at
!    [Fortran-lang Discourse](https://fortran-lang.discourse.group/t/cardanos-solution-of-the-cubic-equation/111/5)
!
!### History
!  * Jacob Williams, 9/14/2022 : created this routine.
!
!@note Works only for single and double precision.

    subroutine polyroots(n, p, wr, wi, info)

    implicit none

    integer,intent(in)   :: n      !! polynomial degree
    real(wp),intent(in)  :: p(n+1) !! polynomial coefficients array, in order of decreasing powers
    real(wp),intent(out) :: wr(n)  !! real parts of roots
    real(wp),intent(out) :: wi(n)  !! imaginary parts of roots
    integer,intent(out)  :: info   !! output from the lapack solver.
                                   !! if `info=0` the routine converged.
                                   !! if `info=-999`, then the leading coefficient is zero.

    integer :: i !! counter
    real(wp),allocatable,dimension(:,:) :: a !! companion matrix
    real(wp),allocatable,dimension(:) :: work !! work array
    real(wp),dimension(1) :: vl, vr !! not used here

#ifdef REAL32
    interface
        subroutine sgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
            implicit none
            character :: jobvl, jobvr
            integer :: info, lda, ldvl, ldvr, lwork, n
            real :: a( lda, * ), vl( ldvl, * ), vr( ldvr, * ), wi( * ), work( * ), wr( * )
        end subroutine sgeev
    end interface
#elif REAL128
    ! do not have a quad solver in lapack
#else
    interface
        subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
            implicit none
            character :: jobvl, jobvr
            integer :: info, lda, ldvl, ldvr, lwork, n
            double precision :: a( lda, * ), vl( ldvl, * ), vr( ldvr, * ), wi( * ), work( * ), wr( * )
        end subroutine dgeev
    end interface
#endif

    ! error check:
    if ( abs(p(1))==0.0_wp ) then
        info = -999
        return
    endif

    ! allocate the work array:
    allocate(work(3*n))

    ! create the companion matrix
    allocate(a(n,n))
    a = 0.0_wp
    do i = 1, n-1
        a(i,i+1) = 1.0_wp
    end do
    do i = n,1,-1
        a(n,n-i+1) = -p(i+1) / p(1)
    end do

    ! call the lapack solver:
#ifdef REAL32
    ! single precision
    call sgeev('N','N',n,a,n,wr,wi,vl,1,vr,1,work,3*n,info)
#elif REAL128
    error stop 'do not have a quad solver in lapack'
#else
    ! by default, use double precision:
    call dgeev('N','N',n,a,n,wr,wi,vl,1,vr,1,work,3*n,info)
#endif

    end subroutine polyroots
!*****************************************************************************************

!*****************************************************************************************
!>
!  This routine computes all zeros of a polynomial of degree NDEG
!  with complex coefficients by computing the eigenvalues of the
!  companion matrix.
!
!### REVISION HISTORY  (YYMMDD)
!  * 791201  DATE WRITTEN. Vandevender, W. H., (SNLA)
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 890531  REVISION DATE from Version 3.2
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!  * 900326  Removed duplicate information from DESCRIPTION section. (WRB)
!  * 911010  Code reworked and simplified.  (RWC and WRB)
!  * Jacob Williams, 9/14/2022 : modernized this code

    subroutine cpqr79(ndeg,coeff,root,ierr)
        implicit none

        integer,intent(in) :: ndeg !! degree of polynomial
        complex(wp),intent(in) :: coeff(*) !! `(NDEG+1)` coefficients in descending order.  i.e.,
                                           !! `P(Z)= COEFF(1)*(Z**NDEG) + COEFF(NDEG)*Z + COEFF(NDEG+1)`
        complex(wp),intent(out) :: root(*) !! `(NDEG)` vector of roots
        integer,intent(out) :: ierr !! Output Error Code.
                                    !!
                                    !!### Normal Code:
                                    !!
                                    !!  * 0 -- means the roots were computed.
                                    !!
                                    !!### Abnormal Codes:
                                    !!
                                    !!  * 1 --  more than 30 QR iterations on some eigenvalue of the companion matrix
                                    !!  * 2 --  COEFF(1)=0.0
                                    !!  * 3 --  NDEG is invalid (less than or equal to 0)

        complex(wp) :: scale , c
        integer :: k , khr , khi , kwr , kwi , kad , kj, km1
        real(wp),dimension(:),allocatable :: work !! work array of dimension `2*NDEG*(NDEG+1)`

        ierr = 0
        if ( abs(coeff(1))==0.0_wp ) then
           ierr = 2
           write(*,*) 'leading coefficient is zero.'
           return
        endif

        if ( ndeg<=0 ) then
           ierr = 3
           write(*,*) 'degree invalid.'
           return
        endif

        if ( ndeg==1 ) then
           root(1) = -coeff(2)/coeff(1)
           return
        endif

        ! allocate work array:
        allocate(work(2*NDEG*(NDEG+1)))

        scale = 1.0_wp/coeff(1)
        khr = 1
        khi = khr + ndeg*ndeg
        kwr = khi + khi - khr
        kwi = kwr + ndeg

        do k = 1 , kwr
           work(k) = 0.0_wp
        enddo

        do k = 1 , ndeg
           kad = (k-1)*ndeg + 1
           c = scale*coeff(k+1)
           work(kad) = -real(c,wp)
           kj = khi + kad - 1
           work(kj) = -aimag(c)
           if ( k/=ndeg ) work(kad+k) = 1.0_wp
        enddo

        call comqr(ndeg,ndeg,1,ndeg,work(khr),work(khi),work(kwr),work(kwi),ierr)

        if ( ierr/=0 ) then
           ierr = 1
           write(*,*) 'no convergence in 30 qr iterations.'
           return
        endif

        do k = 1 , ndeg
           km1 = k - 1
           root(k) = cmplx(work(kwr+km1),work(kwi+km1), wp)
        enddo

        end subroutine cpqr79
!*****************************************************************************************

!*****************************************************************************************
!>
!  this subroutine finds the eigenvalues of a complex
!  upper hessenberg matrix by the qr method.
!
!### Notes
!  * calls [[cdiv]] for complex division.
!  * calls [[csroot]] for complex square root.
!  * calls [[pythag]] for dsqrt(a*a + b*b) .
!
!### Reference
!  * this subroutine is a translation of a unitary analogue of the
!    algol procedure  comlr, num. math. 12, 369-376(1968) by martin
!    and wilkinson.
!    handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
!    the unitary analogue substitutes the qr algorithm of francis
!    (comp. jour. 4, 332-345(1962)) for the lr algorithm.
!
!### History
!  * this version dated august 1983.
!    questions and comments should be directed to burton s. garbow,
!    mathematics and computer science div, argonne national laboratory
!  * Jacob Williams, 9/14/2022 : modernized this code

    subroutine comqr(nm,n,low,igh,hr,hi,wr,wi,ierr)
    implicit none

    integer,intent(in) :: nm !! row dimension of two-dimensional array parameters
    integer,intent(in) :: n !! the order of the matrix
    integer,intent(in) :: low !! integer determined by the balancing
                              !! subroutine  cbal.  if  cbal  has not been used,
                              !! set low=1
    integer,intent(in) :: igh !! integer determined by the balancing
                              !! subroutine  cbal.  if  cbal  has not been used,
                              !! igh=n.
    real(wp),intent(inout) :: hr(nm,n) !! On input: hr and hi contain the real and imaginary parts,
                                       !! respectively, of the complex upper hessenberg matrix.
                                       !! their lower triangles below the subdiagonal contain
                                       !! information about the unitary transformations used in
                                       !! the reduction by  corth, if performed.
                                       !!
                                       !! On Output: the upper hessenberg portions of hr and hi have been
                                       !! destroyed.  therefore, they must be saved before
                                       !! calling  comqr  if subsequent calculation of
                                       !! eigenvectors is to be performed.
    real(wp),intent(inout) :: hi(nm,n) !! See `hr` description
    real(wp),intent(out) :: wr(n) !! the real parts of the eigenvalues.  if an error
                                  !! exit is made, the eigenvalues should be correct
                                  !! for indices `ierr+1,...,n`.
    real(wp),intent(out) :: wi(n) !! the imaginary parts of the eigenvalues.  if an error
                                  !! exit is made, the eigenvalues should be correct
                                  !! for indices `ierr+1,...,n`.
    integer,intent(out) :: ierr !! is set to:
                                !!
                                !!  * 0 -- for normal return
                                !!  * j -- if the limit of 30*n iterations is exhausted
                                !!    while the j-th eigenvalue is being sought.

    integer :: i , j , l , en , ll , itn , its , lp1 , enm1
    real(wp) :: si , sr , ti , tr , xi , xr , yi , yr , zzi , &
                zzr , norm , tst1 , tst2

    ierr = 0
    if ( low/=igh ) then
        ! create real subdiagonal elements
        l = low + 1
        do i = l , igh
            ll = min(i+1,igh)
            if ( hi(i,i-1)/=0.0_wp ) then
                norm = pythag(hr(i,i-1),hi(i,i-1))
                yr = hr(i,i-1)/norm
                yi = hi(i,i-1)/norm
                hr(i,i-1) = norm
                hi(i,i-1) = 0.0_wp
                do j = i , igh
                si = yr*hi(i,j) - yi*hr(i,j)
                hr(i,j) = yr*hr(i,j) + yi*hi(i,j)
                hi(i,j) = si
                enddo
                do j = low , ll
                si = yr*hi(j,i) + yi*hr(j,i)
                hr(j,i) = yr*hr(j,i) - yi*hi(j,i)
                hi(j,i) = si
                enddo
            endif
        enddo
    endif
    ! store roots isolated by cbal
    do i = 1 , n
        if ( i<low .or. i>igh ) then
            wr(i) = hr(i,i)
            wi(i) = hi(i,i)
        endif
    enddo

    en = igh
    tr = 0.0_wp
    ti = 0.0_wp
    itn = 30*n

    ! search for next eigenvalue
100 if ( en<low ) return
    its = 0
    enm1 = en - 1

    ! look for single small sub-diagonal element
    ! for l=en step -1 until low d0 --
200 do ll = low , en
        l = en + low - ll
        if ( l==low ) exit
        tst1 = abs(hr(l-1,l-1)) + abs(hi(l-1,l-1)) + abs(hr(l,l)) + abs(hi(l,l))
        tst2 = tst1 + abs(hr(l,l-1))
        if ( tst2==tst1 ) exit
    enddo

    ! form shift
    if ( l==en ) then
        ! a root found
        wr(en) = hr(en,en) + tr
        wi(en) = hi(en,en) + ti
        en = enm1
        goto 100
    elseif ( itn==0 ) then
        ! set error -- all eigenvalues have not
        !              converged after 30*n iterations
        ierr = en
    else
        if ( its==10 .or. its==20 ) then
            ! form exceptional shift
            sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
            si = 0.0_wp
        else
            sr = hr(en,en)
            si = hi(en,en)
            xr = hr(enm1,en)*hr(en,enm1)
            xi = hi(enm1,en)*hr(en,enm1)
            if ( xr/=0.0_wp .or. xi/=0.0_wp ) then
                yr = (hr(enm1,enm1)-sr)/2.0_wp
                yi = (hi(enm1,enm1)-si)/2.0_wp
                call csroot(yr**2-yi**2+xr,2.0_wp*yr*yi+xi,zzr,zzi)
                if ( yr*zzr+yi*zzi<0.0_wp ) then
                    zzr = -zzr
                    zzi = -zzi
                endif
                call cdiv(xr,xi,yr+zzr,yi+zzi,xr,xi)
                sr = sr - xr
                si = si - xi
            endif
        endif

        do i = low , en
            hr(i,i) = hr(i,i) - sr
            hi(i,i) = hi(i,i) - si
        enddo

        tr = tr + sr
        ti = ti + si
        its = its + 1
        itn = itn - 1
        ! reduce to triangle (rows)
        lp1 = l + 1

        do i = lp1 , en
            sr = hr(i,i-1)
            hr(i,i-1) = 0.0_wp
            norm = pythag(pythag(hr(i-1,i-1),hi(i-1,i-1)),sr)
            xr = hr(i-1,i-1)/norm
            wr(i-1) = xr
            xi = hi(i-1,i-1)/norm
            wi(i-1) = xi
            hr(i-1,i-1) = norm
            hi(i-1,i-1) = 0.0_wp
            hi(i,i-1) = sr/norm

            do j = i , en
                yr = hr(i-1,j)
                yi = hi(i-1,j)
                zzr = hr(i,j)
                zzi = hi(i,j)
                hr(i-1,j) = xr*yr + xi*yi + hi(i,i-1)*zzr
                hi(i-1,j) = xr*yi - xi*yr + hi(i,i-1)*zzi
                hr(i,j) = xr*zzr - xi*zzi - hi(i,i-1)*yr
                hi(i,j) = xr*zzi + xi*zzr - hi(i,i-1)*yi
            enddo

        enddo

        si = hi(en,en)
        if ( si/=0.0_wp ) then
            norm = pythag(hr(en,en),si)
            sr = hr(en,en)/norm
            si = si/norm
            hr(en,en) = norm
            hi(en,en) = 0.0_wp
        endif
        ! inverse operation (columns)
        do j = lp1 , en
            xr = wr(j-1)
            xi = wi(j-1)

            do i = l , j
                yr = hr(i,j-1)
                yi = 0.0_wp
                zzr = hr(i,j)
                zzi = hi(i,j)
                if ( i/=j ) then
                    yi = hi(i,j-1)
                    hi(i,j-1) = xr*yi + xi*yr + hi(j,j-1)*zzi
                endif
                hr(i,j-1) = xr*yr - xi*yi + hi(j,j-1)*zzr
                hr(i,j) = xr*zzr + xi*zzi - hi(j,j-1)*yr
                hi(i,j) = xr*zzi - xi*zzr - hi(j,j-1)*yi
            enddo

        enddo

        if ( si/=0.0_wp ) then
            do i = l , en
                yr = hr(i,en)
                yi = hi(i,en)
                hr(i,en) = sr*yr - si*yi
                hi(i,en) = sr*yi + si*yr
            enddo
        endif
        goto 200
    endif

    end subroutine comqr
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the complex square root of a complex number.
!
!  `(YR,YI) = complex sqrt(XR,XI)`
!
!### REVISION HISTORY  (YYMMDD)
!  * 811101  DATE WRITTEN
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900402  Added TYPE section.  (WRB)
!  * Jacob Williams, 9/14/2022 : modernized this code

    pure subroutine csroot(xr,xi,yr,yi)
    implicit none

    real(wp),intent(in) :: xr , xi
    real(wp),intent(out) :: yr , yi

    real(wp) :: s , tr , ti

    ! branch chosen so that yr >= 0.0 and sign(yi) == sign(xi)
    tr = xr
    ti = xi
    s = sqrt(0.5_wp*(pythag(tr,ti)+abs(tr)))
    if ( tr>=0.0_wp ) yr = s
    if ( ti<0.0_wp ) s = -s
    if ( tr<=0.0_wp ) yi = s
    if ( tr<0.0_wp ) yr = 0.5_wp*(ti/yi)
    if ( tr>0.0_wp ) yi = 0.5_wp*(ti/yr)

    end subroutine csroot
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the complex quotient of two complex numbers.
!
!  Complex division, `(CR,CI) = (AR,AI)/(BR,BI)`
!
!### REVISION HISTORY  (YYMMDD)
!  * 811101  DATE WRITTEN
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900402  Added TYPE section.  (WRB)
!  * Jacob Williams, 9/14/2022 : modernized this code

    pure subroutine cdiv(ar,ai,br,bi,cr,ci)
    implicit none

    real(wp),intent(in) :: ar , ai , br , bi
    real(wp),intent(out) :: cr , ci

    real(wp) :: s , ars , ais , brs , bis

    s = abs(br) + abs(bi)
    ars = ar/s
    ais = ai/s
    brs = br/s
    bis = bi/s
    s = brs**2 + bis**2
    cr = (ars*brs+ais*bis)/s
    ci = (ais*brs-ars*bis)/s

    end subroutine cdiv
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the complex square root of a complex number without
!  destructive overflow or underflow.
!
!  Finds `sqrt(A**2+B**2)` without overflow or destructive underflow
!
!### REVISION HISTORY  (YYMMDD)
!  * 811101  DATE WRITTEN
!  * 890531  Changed all specific intrinsics to generic.  (WRB)
!  * 891214  Prologue converted to Version 4.0 format.  (BAB)
!  * 900402  Added TYPE section.  (WRB)
!  * Jacob Williams, 9/14/2022 : modernized this code

    pure real(wp) function pythag(a,b)
    implicit none

    real(wp),intent(in) :: a , b

    real(wp) :: p , q , r , s , t

    p = max(abs(a),abs(b))
    q = min(abs(a),abs(b))

    if ( q==0.0_wp ) then
        pythag = p
    else
        do
            r = (q/p)**2
            t = 4.0_wp + r
            if ( t==4.0_wp ) then
                pythag = p
                exit
            else
                s = r/t
                p = p + 2.0_wp*p*s
                q = q*s
            endif
        end do
    endif

    end function pythag
!*****************************************************************************************

!*****************************************************************************************
!>
!  "Polish" a root using a complex Newton Raphson method.
!  This routine will perform a Newton iteration and update the roots only if they improve,
!  otherwise, they are left as is.
!
!### History
!  * Jacob Williams, 9/15/2023, created this routine.

    subroutine newton_root_polish(n, p, zr, zi, ftol, ztol, maxiter, istat)

    implicit none

    integer, intent(in)     :: n         !! degree of polynomial
    real(wp), intent(in)    :: p(n+1)    !! vector of coefficients in order of decreasing powers
    real(wp), intent(inout) :: zr        !! output vector of real parts of the zeros
    real(wp), intent(inout) :: zi        !! output vector of imaginary parts of the zeros
    real(wp), intent(in)    :: ftol      !! convergence tolerance for the root
    real(wp), intent(in)    :: ztol      !! convergence tolerance for `x`
    integer, intent(in)     :: maxiter   !! maximum number of iterations
    integer, intent(out)    :: istat     !! status flag:
                                         !!
                                         !! * 0  = converged in `f`
                                         !! * 1  = converged in `x`
                                         !! * -1 = singular
                                         !! * -2 = max iterations reached

    complex(wp) :: z, f, z_prev, z_best, f_best, dfdx
    integer :: i !! counter

    real(wp),parameter :: alpha = 1.0_wp !! newton step size

    ! first evaluate initial point:
    z = cmplx(zr, zi, wp)
    call func(z,f,dfdx)

    ! initialize:
    istat = 0
    z_prev = z
    z_best = z
    f_best = f

    main : do i = 1, maxiter

        if (i>1) call func(z,f,dfdx)
        if (abs(f)<abs(f_best)) then
            ! best so far
            zr = real(z, wp)
            zi = aimag(z)
            z_best = z
            f_best = f
        end if
        if (abs(f)<=ftol) exit main

        if (i == maxiter) then ! max iterations reached
            istat = -2
            exit main
        end if

        if (dfdx==0.0_wp) then  ! can't proceed
            istat = -1
            exit main
        end if

        ! Newton correction for next step:
        z = z - alpha * ( f / dfdx )

        if (abs(z-z_prev)<=ztol) then
            ! convergence in x. try this point and see if there is any improvement
            istat = 1
            call func(z,f,dfdx)
            if (abs(f)<abs(f_best)) then ! best so far
                zr = real(z, wp)
                zi = aimag(z)
            end if
            exit main
        end if
        z_prev = z

    end do main

    contains

        subroutine func(x,f,dfdx)

            !! evaluate the polynomial:
            !! `f = p(1)*x**n + p(2)*x**(n-1) + ... + p(n)*x + p(n+1)`
            !! and its derivative using Horner's Rule.
            !!
            !! See: "Roundoff in polynomial evaluation", W. Kahan, 1986
            !! https://rosettacode.org/wiki/Horner%27s_rule_for_polynomial_evaluation#Fortran

            complex(wp),intent(in)  :: x
            complex(wp),intent(out) :: f    !! function value at x
            complex(wp),intent(out) :: dfdx !! function derivative at x

            integer :: i !! counter

            f = p(1)
            dfdx = 0.0_wp
            do i = 2, n+1
                dfdx = dfdx * x + f
                f = f * x + p(i)
            end do

        end subroutine func

    end subroutine newton_root_polish
!*****************************************************************************************

!*****************************************************************************************
    end module polyroots_module
!*****************************************************************************************