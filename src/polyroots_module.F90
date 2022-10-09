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
    use ieee_arithmetic

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
    real(wp), parameter :: deg2rad = pi/180.0_wp

    ! general polynomial routines:
    public :: polyroots
    public :: cpolyroots
    public :: rpoly
    public :: cpoly
    public :: cpzero
    public :: rpzero
    public :: rpqr79
    public :: cpqr79
    public :: qr_algeq_solver
    public :: cmplx_roots_gen
    public :: polzeros
    public :: fpml
    public :: dpolz
    public :: cpolz

    ! special polynomial routines:
    public :: dcbcrt
    public :: dqdcrt
    public :: rroots_chebyshev_cubic

    ! utility routines:
    public :: newton_root_polish
    public :: sort_roots

    interface newton_root_polish
        module procedure :: newton_root_polish_real, &
                            newton_root_polish_complex
    end interface

contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finds the zeros of a general real polynomial using the Jenkins & Traub algorithm.
!
!### History
!  * M. A. Jenkins, "[Algorithm 493: Zeros of a Real Polynomial](https://dl.acm.org/doi/10.1145/355637.355643)",
!    ACM Transactions on Mathematical SoftwareVolume 1, Issue 2, June 1975, pp 178-189
!  * code converted using to_f90 by alan miller, 2003-06-02
!  * Jacob Williams, 9/13/2022 : modernized this code

subroutine rpoly(op, degree, zeror, zeroi, istat)

    implicit none

    integer, intent(in) :: degree !! degree of polynomial
    real(wp), dimension(degree+1), intent(in) :: op !! vector of coefficients in order of decreasing powers
    real(wp), dimension(degree), intent(out) :: zeror !! output vector of real parts of the zeros
    real(wp), dimension(degree), intent(out) :: zeroi !! output vector of imaginary parts of the zeros
    integer, intent(out) :: istat !! status output:
                                  !!
                                  !! * `0` : success
                                  !! * `-1` : leading coefficient is zero
                                  !! * `-2` : no roots found
                                  !! * `>0` : the number of zeros found

    real(wp), dimension(:), allocatable :: p, qp, k, qk, svk, temp, pt
    real(wp) :: sr, si, u, v, a, b, c, d, a1, a3, &
                a7, e, f, g, h, szr, szi, lzr, lzi, &
                t, aa, bb, cc, factor, mx, mn, xx, yy, &
                xxx, x, sc, bnd, xm, ff, df, dx
    integer :: cnt, nz, i, j, jj, l, nm1, n, nn
    logical :: zerok

    real(wp), parameter :: cosr = cos(94.0_wp*deg2rad)
    real(wp), parameter :: sinr = sin(86.0_wp*deg2rad)
    real(wp), parameter :: base = radix(1.0_wp)
    real(wp), parameter :: eta = eps
    real(wp), parameter :: infin = huge(1.0_wp)
    real(wp), parameter :: smalno = tiny(1.0_wp)
    real(wp), parameter :: sqrthalf = sqrt(0.5_wp)
    real(wp), parameter :: are = eta !! unit error in +
    real(wp), parameter :: mre = eta !! unit error in *
    real(wp), parameter :: lo = smalno/eta

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

    main: do

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
        scale: block
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
    if (istat == 0) istat = -2  ! if not roots found

contains

    subroutine fxshfr(l2, nz)

      !! computes up to  l2  fixed shift k-polynomials, testing for convergence in
      !! the linear or quadratic case.  initiates one of the variable shift
      !! iterations and returns with the number of zeros found.

        integer, intent(in) :: l2 !! limit of fixed shift steps
        integer, intent(out) :: nz !! number of zeros found

        real(wp) :: svu, svv, ui, vi, s, betas, betav, oss, ovv, &
                    ss, vv, ts, tv, ots, otv, tvv, tss
        integer :: type, j, iflag
        logical :: vpass, spass, vtry, stry, skip

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

                        try: block

                            if (.not. skip) then
                                call quadit(ui, vi, nz)
                                if (nz > 0) return

                                ! quadratic iteration has failed. flag that it has
                                ! been tried and decrease the convergence criterion.
                                vtry = .true.
                                betav = betav*0.25_wp

                                ! try linear iteration if it has not been tried and
                                ! the s sequence is converging
                                if (stry .or. (.not. spass)) exit try
                                k(1:n) = svk(1:n)
                            end if
                            skip = .false.

                            call realit(s, nz, iflag)
                            if (nz > 0) return

                            ! linear iteration has failed.  flag that it has been
                            ! tried and decrease the convergence criterion
                            stry = .true.
                            betas = betas*0.25_wp
                            if (iflag /= 0) then
                                ! if linear iteration signals an almost double real
                                ! zero attempt quadratic interation
                                ui = -(s + s)
                                vi = s*s
                                cycle
                            end if

                        end block try

                        ! restore variables
                        u = svu
                        v = svv
                        k(1:n) = svk(1:n)

                        ! try quadratic iteration if it has not been tried
                        ! and the v sequence is converging
                        if (.not. (vpass .and. (.not. vtry))) exit

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

        real(wp), intent(in) :: uu !! coefficients of starting quadratic
        real(wp), intent(in) :: vv !! coefficients of starting quadratic
        integer, intent(out) :: nz !! number of zero found

        real(wp) :: ui, vi, mp, omp, ee, relstp, t, zm
        integer :: type, i, j
        logical :: tried

        nz = 0
        tried = .false.
        u = uu
        v = vv
        j = 0

        ! main loop
        main: do
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
        integer, intent(out) :: nz !! number of zero found
        integer, intent(out) :: iflag !! flag to indicate a pair of zeros near real axis.

        real(wp) :: pv, kv, t, s, ms, mp, omp, ee
        integer :: i, j

        nz = 0
        s = sss
        iflag = 0
        j = 0

        ! main loop
        main: do
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
                if (abs(t) <= 0.001_wp*abs(s - t) .and. mp > omp) then
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

        integer, intent(out) :: type !! integer variable set here indicating how the
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
        else
            type = 1
            ! type=1 indicates that all formulas are divided by c
            e = a/c
            f = d/c
            g = u*e
            h = v*b
            a3 = a*e + (h/c + g)*b
            a1 = b - a*(d/c)
            a7 = a + g*d + h*f
        end if

    end subroutine calcsc

    subroutine nextk(type)

         !! computes the next k polynomials using scalars computed in calcsc.

        integer, intent(in) :: type

        real(wp) :: temp
        integer :: i

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

        else
            ! use unscaled form of the recurrence if type is 3
            k(1) = 0.0_wp
            k(2) = 0.0_wp
            do i = 3, n
                k(i) = qk(i - 2)
            end do
        end if

    end subroutine nextk

    subroutine newest(type, uu, vv)

        ! compute new estimates of the quadratic coefficients
        ! using the scalars computed in calcsc.

        integer, intent(in) :: type
        real(wp), intent(out) :: uu
        real(wp), intent(out) :: vv

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
                vv = v*(1.0_wp + c4/temp)
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

        integer, intent(in) :: nn
        real(wp), intent(in) :: u, v, p(nn)
        real(wp), intent(out) :: q(nn), a, b

        real(wp) :: c
        integer :: i

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

        real(wp), intent(in) :: a, b1, c
        real(wp), intent(out) :: sr, si, lr, li

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

    real(wp), dimension(4), intent(in) :: a     !! coefficients
    real(wp), dimension(3), intent(out) :: zr   !! real components of roots
    real(wp), dimension(3), intent(out) :: zi   !! imaginary components of roots

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

    real(wp), dimension(3), intent(in) :: a     !! coefficients
    real(wp), dimension(2), intent(out) :: zr   !! real components of roots
    real(wp), dimension(2), intent(out) :: zi   !! imaginary components of roots

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

subroutine qr_algeq_solver(n, c, zr, zi, istatus, detil)

    implicit none

    integer, intent(in) :: n !! degree of the monic polynomial.
    real(wp), intent(in) :: c(n + 1) !! coefficients of the polynomial. in order of decreasing powers.
    real(wp), intent(out) :: zr(n) !! real part of output roots
    real(wp), intent(out) :: zi(n) !! imaginary part of output roots
    integer, intent(out) :: istatus !! return code:
                                    !!
                                    !! * -1 : degree <= 0
                                    !! * -2 : leading coefficient `c(1)` is zero
                                    !! * 0 : success
                                    !! * otherwise, the return code from `hqr_eigen_hessenberg`
    real(wp), intent(out), optional :: detil !! accuracy hint.

    real(wp), allocatable :: a(:, :) !! work matrix
    integer, allocatable :: cnt(:) !! work area for counting the qr-iterations
    integer :: i, iter
    real(wp) :: afnorm

    ! check for invalid arguments
    if (n <= 0) then
        istatus = -1
        return
    end if
    if (c(1) == 0.0_wp) then
        ! leading coefficient is zero.
        istatus = -2
        return
    end if

    allocate (a(n, n))
    allocate (cnt(n))

    ! build the companion matrix a.
    call build_companion(n, a, c)

    ! balancing the a itself.
    call balance_companion(n, a)

    ! qr iterations from a.
    call hqr_eigen_hessenberg(n, a, zr, zi, cnt, istatus)
    if (istatus /= 0) then
        write (*, '(A,1X,I4)') 'abnormal return from hqr_eigen_hessenberg. istatus=', istatus
        if (istatus == 1) write (*, '(A)') 'matrix is completely zero.'
        if (istatus == 2) write (*, '(A)') 'qr iteration did not converge.'
        if (istatus > 3) write (*, '(A)') 'arguments violate the condition.'
        return
    end if

    if (present(detil)) then

        ! compute the frobenius norm of the balanced companion matrix a.
        afnorm = frobenius_norm_companion(n, a)

        ! count the total qr iteration.
        iter = 0
        do i = 1, n
            if (cnt(i) > 0) iter = iter + cnt(i)
        end do

        ! calculate the accuracy hint.
        detil = eps*n*iter*afnorm

    end if

contains

    subroutine build_companion(n, a, c)

        !!  build the companion matrix of the polynomial.
        !!  (this was modified to allow for non-monic polynomials)

        implicit none

        integer, intent(in) :: n
        real(wp), intent(out) :: a(n, n)
        real(wp), intent(in) :: c(n + 1) !! coefficients in order of decreasing powers

        integer :: i !! counter

        ! create the companion matrix
        a = 0.0_wp
        do i = 1, n - 1
            a(i + 1, i) = 1.0_wp
        end do
        do i = n, 1, -1
            a(n - i + 1, n) = -c(i + 1)/c(1)
        end do

    end subroutine build_companion

    subroutine balance_companion(n, a)

        !!  blancing the unsymmetric matrix `a`.
        !!
        !!  this fortran code is based on the algol code "balance" from paper:
        !!   "balancing a matrix for calculation of eigenvalues and eigenvectors"
        !!   by b.n.parlett and c.reinsch, numer. math. 13, 293-304(1969).
        !!
        !!  note: the only non-zero elements of the companion matrix are touched.

        implicit none

        integer, intent(in) :: n
        real(wp), intent(inout) :: a(n, n)

        integer, parameter :: b = radix(1.0_wp) !! base of the floating point representation on the machine
        integer, parameter :: b2 = b**2

        integer :: i, j
        real(wp) :: c, f, g, r, s
        logical :: noconv

        if (n <= 1) return ! do nothing

        ! iteration:
        main: do
            noconv = .false.
            do i = 1, n
                ! touch only non-zero elements of companion.
                if (i /= n) then
                    c = abs(a(i + 1, i))
                else
                    c = 0.0_wp
                    do j = 1, n - 1
                        c = c + abs(a(j, n))
                    end do
                end if
                if (i == 1) then
                    r = abs(a(1, n))
                elseif (i /= n) then
                    r = abs(a(i, i - 1)) + abs(a(i, n))
                else
                    r = abs(a(i, i - 1))
                end if

                if (c /= 0.0_wp .and. r /= 0.0_wp) then

                    g = r/b
                    f = 1.0_wp
                    s = c + r
                    do
                        if (c >= g) exit
                        f = f*b
                        c = c*b2
                    end do
                    g = r*b
                    do
                        if (c < g) exit
                        f = f/b
                        c = c/b2
                    end do
                    if ((c + r)/f < 0.95_wp*s) then
                        g = 1.0_wp/f
                        noconv = .true.
                        ! touch only non-zero elements of companion.
                        if (i == 1) then
                            a(1, n) = a(1, n)*g
                        else
                            a(i, i - 1) = a(i, i - 1)*g
                            a(i, n) = a(i, n)*g
                        end if
                        if (i /= n) then
                            a(i + 1, i) = a(i + 1, i)*f
                        else
                            do j = 1, n
                                a(j, i) = a(j, i)*f
                            end do
                        end if
                    end if
                end if
            end do
            if (noconv) cycle main
            exit main
        end do main

    end subroutine balance_companion

    function frobenius_norm_companion(n, a) result(afnorm)

        !!  calculate the frobenius norm of the companion-like matrix.
        !!  note: the only non-zero elements of the companion matrix are touched.

        implicit none

        integer, intent(in) :: n
        real(wp), intent(in) :: a(n, n)
        real(wp) :: afnorm

        integer :: i, j
        real(wp) :: sum

        sum = 0.0_wp
        do j = 1, n - 1
            sum = sum + a(j + 1, j)**2
        end do
        do i = 1, n
            sum = sum + a(i, n)**2
        end do
        afnorm = sqrt(sum)

    end function frobenius_norm_companion

    subroutine hqr_eigen_hessenberg(n0, h, wr, wi, cnt, istatus)

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
        !!          `maxiter` iterations.

        implicit none

        integer, intent(in) :: n0
        real(wp), intent(inout) :: h(n0, n0)
        real(wp), intent(out) :: wr(n0)
        real(wp), intent(out) :: wi(n0)
        integer, intent(inout) :: cnt(n0)
        integer, intent(out) :: istatus

        integer :: i, j, k, l, m, na, its, n
        real(wp) :: p, q, r, s, t, w, x, y, z
        logical :: notlast, found

#if REAL128
        integer, parameter :: maxiter = 100 !! max iterations. It seems we need more than 30
                                            !! for quad precision (see test case 11)
#else
        integer, parameter :: maxiter = 30  !! max iterations
#endif

        ! note: n is changing in this subroutine.
        n = n0
        istatus = 0
        t = 0.0_wp

        main: do

            if (n == 0) return

            its = 0
            na = n - 1

            do

                ! look for single small sub-diagonal element
                found = .false.
                do l = n, 2, -1
                    if (abs(h(l, l - 1)) <= eps*(abs(h(l - 1, l - 1)) + abs(h(l, l)))) then
                        found = .true.
                        exit
                    end if
                end do
                if (.not. found) l = 1

                x = h(n, n)
                if (l == n) then
                    ! one root found
                    wr(n) = x + t
                    wi(n) = 0.0_wp
                    cnt(n) = its
                    n = na
                    cycle main
                else
                    y = h(na, na)
                    w = h(n, na)*h(na, n)
                    if (l == na) then
                        ! comment: two roots found
                        p = (y - x)/2
                        q = p**2 + w
                        y = sqrt(abs(q))
                        cnt(n) = -its
                        cnt(na) = its
                        x = x + t
                        if (q > 0.0_wp) then
                            ! real pair
                            if (p < 0.0_wp) y = -y
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
                        end if
                        n = n - 2
                        cycle main
                    else
                        if (its == maxiter) then ! 30 for the original double precision code
                            istatus = 1
                            return
                        end if
                        if (its == 10 .or. its == 20) then
                            ! form exceptional shift
                            t = t + x
                            do i = 1, n
                                h(i, i) = h(i, i) - x
                            end do
                            s = abs(h(n, na)) + abs(h(na, n - 2))
                            y = 0.75_wp*s
                            x = y
                            w = -0.4375_wp*s**2
                        end if
                        its = its + 1
                        ! look for two consecutive small sub-diagonal elements
                        do m = n - 2, l, -1
                            z = h(m, m)
                            r = x - z
                            s = y - z
                            p = (r*s - w)/h(m + 1, m) + h(m, m + 1)
                            q = h(m + 1, m + 1) - z - r - s
                            r = h(m + 2, m + 1)
                            s = abs(p) + abs(q) + abs(r)
                            p = p/s
                            q = q/s
                            r = r/s
                            if (m == l) exit
                            if (abs(h(m, m - 1))*(abs(q) + abs(r)) <= eps*abs(p) &
                                *(abs(h(m - 1, m - 1)) + abs(z) + abs(h(m + 1, m + 1)))) exit
                        end do

                        do i = m + 2, n
                            h(i, i - 2) = 0.0_wp
                        end do
                        do i = m + 3, n
                            h(i, i - 3) = 0.0_wp
                        end do
                        ! double qr step involving rows l to n and columns m to n
                        do k = m, na
                            notlast = (k /= na)
                            if (k /= m) then
                                p = h(k, k - 1)
                                q = h(k + 1, k - 1)
                                if (notlast) then
                                    r = h(k + 2, k - 1)
                                else
                                    r = 0.0_wp
                                end if
                                x = abs(p) + abs(q) + abs(r)
                                if (x == 0.0_wp) cycle
                                p = p/x
                                q = q/x
                                r = r/x
                            end if
                            s = sqrt(p**2 + q**2 + r**2)
                            if (p < 0.0_wp) s = -s
                            if (k /= m) then
                                h(k, k - 1) = -s*x
                            elseif (l /= m) then
                                h(k, k - 1) = -h(k, k - 1)
                            end if
                            p = p + s
                            x = p/s
                            y = q/s
                            z = r/s
                            q = q/p
                            r = r/p
                            ! row modification
                            do j = k, n
                                p = h(k, j) + q*h(k + 1, j)
                                if (notlast) then
                                    p = p + r*h(k + 2, j)
                                    h(k + 2, j) = h(k + 2, j) - p*z
                                end if
                                h(k + 1, j) = h(k + 1, j) - p*y
                                h(k, j) = h(k, j) - p*x
                            end do
                            if (k + 3 < n) then
                                j = k + 3
                            else
                                j = n
                            end if
                            ! column modification;
                            do i = l, j
                                p = x*h(i, k) + y*h(i, k + 1)
                                if (notlast) then
                                    p = p + z*h(i, k + 2)
                                    h(i, k + 2) = h(i, k + 2) - p*r
                                end if
                                h(i, k + 1) = h(i, k + 1) - p*q
                                h(i, k) = h(i, k) - p
                            end do
                        end do
                        cycle
                    end if
                end if

            end do

        end do main

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

subroutine cpevl(n, m, a, z, c, b, kbd)

    implicit none

    integer, intent(in) :: n !! Degree of the polynomial
    integer, intent(in) :: m !! Number of derivatives to be calculated:
                             !!
                             !!  * `M=0` evaluates only the function
                             !!  * `M=1` evaluates the function and first derivative, etc.
                             !!
                             !! if `M > N+1` function and all `N` derivatives will be calculated.
    complex(wp), intent(in) :: a(*) !! vector containing the `N+1` coefficients of polynomial.
                                    !! `A(I)` = coefficient of `Z**(N+1-I)`
    complex(wp), intent(in) :: z !! point at which the evaluation is to take place
    complex(wp), intent(out) :: c(*) !! Array of `2(M+1)` words: `C(I+1)` contains the complex value of the I-th
                                     !! derivative at `Z, I=0,...,M`
    complex(wp), intent(out) :: b(*) !! Array of `2(M+1)` words: `B(I)` contains the bounds on the real and imaginary parts
                                     !! of `C(I)` if they were requested. only needed if bounds are to be calculated.
                                     !! It is not used otherwise.
    logical, intent(in) :: kbd !! A logical variable, e.g. .TRUE. or .FALSE. which is
                               !! to be set .TRUE. if bounds are to be computed.

    real(wp) :: r, s
    integer :: i, j, mini, np1
    complex(wp) :: ci, cim1, bi, bim1, t, za, q

    za(q) = cmplx(abs(real(q, wp)), abs(aimag(q)), wp)
    np1 = n + 1
    do j = 1, np1
        ci = 0.0_wp
        cim1 = a(j)
        bi = 0.0_wp
        bim1 = 0.0_wp
        mini = min(m + 1, n + 2 - j)
        do i = 1, mini
            if (j /= 1) ci = c(i)
            if (i /= 1) cim1 = c(i - 1)
            c(i) = cim1 + z*ci
            if (kbd) then
                if (j /= 1) bi = b(i)
                if (i /= 1) bim1 = b(i - 1)
                t = bi + (3.0_wp*eps + 4.0_wp*eps*eps)*za(ci)
                r = real(za(z)*cmplx(real(t, wp), -aimag(t), wp), wp)
                s = aimag(za(z)*t)
                b(i) = (1.0_wp + 8.0_wp*eps)*(bim1 + eps*za(cim1) + cmplx(r, s, wp))
                if (j == 1) b(i) = 0.0_wp
            end if
        end do
    end do

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

subroutine cpzero(in, a, r, iflg, s)

    implicit none

    integer, intent(in) :: in !! `N`, the degree of `P(Z)`
    complex(wp), dimension(in+1), intent(in) :: a !! complex vector containing coefficients of `P(Z)`,
                                                  !! `A(I)` = coefficient of `Z**(N+1-i)`
    complex(wp), dimension(in), intent(inout) :: r !! `N` word complex vector. On input: containing initial
                                                   !! estimates for zeros if these are known. On output: Ith zero
    integer, intent(inout) :: iflg !!### On Input:
                                   !!
                                   !! flag to indicate if initial estimates of zeros are input:
                                   !!
                                   !!  * If `IFLG == 0`, no estimates are input.
                                   !!  * If `IFLG /= 0`, the vector `R` contains estimates of the zeros
                                   !!
                                   !! ** WARNING ****** If estimates are input, they must
                                   !!                   be separated, that is, distinct or
                                   !!                   not repeated.
                                   !!### On Output:
                                   !!
                                   !! Error Diagnostics:
                                   !!
                                   !! * If `IFLG == 0` on return, all is well
                                   !! * If `IFLG == 1` on return, `A(1)=0.0` or `N=0` on input
                                   !! * If `IFLG == 2` on return, the program failed to converge
                                   !!   after `25*N` iterations.  Best current estimates of the
                                   !!   zeros are in `R(I)`.  Error bounds are not calculated.
    real(wp), intent(out) :: s(in) !! an `N` word array. `S(I)` = bound for `R(I)`

    integer :: i, imax, j, n, n1, nit, nmax, nr
    real(wp) :: u, v, x
    complex(wp) :: pn, temp
    complex(wp) :: ctmp(1), btmp(1)
    complex(wp), dimension(:), allocatable :: t !! `4(N+1)` word array used for temporary storage

    if (in <= 0 .or. abs(a(1)) == 0.0_wp) then
        iflg = 1
    else
        ! work array:
        allocate(t(4*(in+1)))
        ! check for easily obtained zeros
        n = in
        n1 = n + 1
        if (iflg == 0) then
            do
                n1 = n + 1
                if (n <= 1) then
                    r(1) = -a(2)/a(1)
                    s(1) = 0.0_wp
                    return
                elseif (abs(a(n1)) /= 0.0_wp) then
                    ! if initial estimates for zeros not given, find some
                    temp = -a(2)/(a(1)*n)
                    call cpevl(n, n, a, temp, t, t, .false.)
                    imax = n + 2
                    t(n1) = abs(t(n1))
                    do i = 2, n1
                        t(n + i) = -abs(t(n + 2 - i))
                        if (real(t(n + i), wp) < real(t(imax), wp)) imax = n + i
                    end do
                    x = (-real(t(imax), wp)/real(t(n1), wp))**(1.0_wp/(imax - n1))
                    do
                        x = 2.0_wp*x
                        call cpevl(n, 0, t(n1), cmplx(x, 0.0_wp, wp), ctmp, btmp, .false.)
                        pn = ctmp(1)
                        if (real(pn, wp) >= 0.0_wp) exit
                    end do
                    u = 0.5_wp*x
                    v = x
                    do
                        x = 0.5_wp*(u + v)
                        call cpevl(n, 0, t(n1), cmplx(x, 0.0_wp, wp), ctmp, btmp, .false.)
                        pn = ctmp(1)
                        if (real(pn, wp) > 0.0_wp) v = x
                        if (real(pn, wp) <= 0.0_wp) u = x
                        if ((v - u) <= 0.001_wp*(1.0_wp + v)) exit
                    end do
                    do i = 1, n
                        u = (pi/n)*(2*i - 1.5_wp)
                        r(i) = max(x, 0.001_wp*abs(temp))*cmplx(cos(u), sin(u), wp) + temp
                    end do
                    exit
                else
                    r(n) = 0.0_wp
                    s(n) = 0.0_wp
                    n = n - 1
                end if
            end do
        end if

        ! main iteration loop starts here
        nr = 0
        nmax = 25*n
        do nit = 1, nmax
            do i = 1, n
                if (nit == 1 .or. abs(t(i)) /= 0.0_wp) then
                    call cpevl(n, 0, a, r(i), ctmp, btmp, .true.)
                    pn = ctmp(1)
                    temp = btmp(1)
                    if (abs(real(pn, wp)) + abs(aimag(pn)) > real(temp, wp) + aimag(temp)) then
                        temp = a(1)
                        do j = 1, n
                            if (j /= i) temp = temp*(r(i) - r(j))
                        end do
                        t(i) = pn/temp
                    else
                        t(i) = 0.0_wp
                        nr = nr + 1
                    end if
                end if
            end do
            do i = 1, n
                r(i) = r(i) - t(i)
            end do
            if (nr == n) then
                ! calculate error bounds for zeros
                do nr = 1, n
                    call cpevl(n, n, a, r(nr), t, t(n + 2), .true.)
                    x = abs(cmplx(abs(real(t(1), wp)), abs(aimag(t(1))), wp) + t(n + 2))
                    s(nr) = 0.0_wp
                    do i = 1, n
                        x = x*real(n1 - i, wp)/i
                        temp = cmplx(max(abs(real(t(i + 1), wp)) - real(t(n1 + i), wp), 0.0_wp), &
                                     max(abs(aimag(t(i + 1))) - aimag(t(n1 + i)), 0.0_wp), wp)
                        s(nr) = max(s(nr), (abs(temp)/x)**(1.0_wp/i))
                    end do
                    s(nr) = 1.0_wp/s(nr)
                end do
                return
            end if
        end do
        iflg = 2  ! error exit
    end if

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

subroutine rpzero(n, a, r, iflg, s)

    implicit none

    integer, intent(in) :: n !! degree of `P(X)`
    real(wp), dimension(n+1), intent(in) :: a !! real vector containing coefficients of `P(X)`,
                                              !! `A(I)` = coefficient of `X**(N+1-I)`
    complex(wp), dimension(n), intent(inout) :: r !! `N` word complex vector. On Input: containing initial estimates for zeros
                                                  !! if these are known. On output: ith zero.
    integer, intent(inout) :: iflg !!### On Input:
                                   !!
                                   !! flag to indicate if initial estimates of zeros are input:
                                   !!
                                   !!  * If `IFLG == 0`, no estimates are input.
                                   !!  * If `IFLG /= 0`, the vector R contains estimates of the zeros
                                   !!
                                   !! ** WARNING ****** If estimates are input, they must
                                   !!                   be separated, that is, distinct or
                                   !!                   not repeated.
                                   !!### On Output:
                                   !!
                                   !! Error Diagnostics:
                                   !!
                                   !! * If `IFLG == 0` on return, all is well
                                   !! * If `IFLG == 1` on return, `A(1)=0.0` or `N=0` on input
                                   !! * If `IFLG == 2` on return, the program failed to converge
                                   !!   after `25*N` iterations.  Best current estimates of the
                                   !!   zeros are in `R(I)`.  Error bounds are not calculated.
    real(wp), dimension(n), intent(out) :: s !! an `N` word array. bound for `R(I)`.

    integer :: i
    complex(wp), dimension(:), allocatable :: p !! complex coefficients

    allocate(p(n+1))

    do i = 1, n + 1
        p(i) = cmplx(a(i), 0.0_wp, wp)
    end do
    call cpzero(n, p, r, iflg, s)

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

subroutine rpqr79(ndeg, coeff, root, ierr)

    implicit none

    integer, intent(in) :: ndeg !! degree of polynomial
    real(wp), dimension(ndeg+1), intent(in) :: coeff !! `NDEG+1` coefficients in descending order.  i.e.,
                                                     !! `P(Z) = COEFF(1)*(Z**NDEG) + COEFF(NDEG)*Z + COEFF(NDEG+1)`
    complex(wp), dimension(ndeg), intent(out) :: root !! `NDEG` vector of roots
    integer, intent(out) :: ierr !! Output Error Code
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
    integer :: k, kh, kwr, kwi, kcol, km1, kwend
    real(wp), dimension(:), allocatable :: work !! work array of dimension `NDEG*(NDEG+2)`

    ierr = 0
    if (abs(coeff(1)) == 0.0_wp) then
        ierr = 2
        write (*, *) 'leading coefficient is zero.'
        return
    end if

    if (ndeg <= 0) then
        ierr = 3
        write (*, *) 'degree invalid.'
        return
    end if

    if (ndeg == 1) then
        root(1) = cmplx(-coeff(2)/coeff(1), 0.0_wp, wp)
        return
    end if

    allocate (work(ndeg*(ndeg + 2))) ! work array

    scale = 1.0_wp/coeff(1)
    kh = 1
    kwr = kh + ndeg*ndeg
    kwi = kwr + ndeg
    kwend = kwi + ndeg - 1

    do k = 1, kwend
        work(k) = 0.0_wp
    end do

    do k = 1, ndeg
        kcol = (k - 1)*ndeg + 1
        work(kcol) = -coeff(k + 1)*scale
        if (k /= ndeg) work(kcol + k) = 1.0_wp
    end do

    call hqr(ndeg, ndeg, 1, ndeg, work(kh), work(kwr), work(kwi), ierr)

    if (ierr /= 0) then
        ierr = 1
        write (*, *) 'no convergence in 30 qr iterations.'
        return
    end if

    do k = 1, ndeg
        km1 = k - 1
        root(k) = cmplx(work(kwr + km1), work(kwi + km1), wp)
    end do

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

subroutine hqr(nm, n, low, igh, h, wr, wi, ierr)

    implicit none

    integer, intent(in) :: nm !! must be set to the row dimension of two-dimensional
                              !! array parameters as declared in the calling program
                              !! dimension statement.
    integer, intent(in) :: n !! order of the matrix
    integer, intent(in) :: low  !! low and igh are integers determined by the balancing
                                !! subroutine  balanc.  if  balanc  has not been used,
                                !! set low=1, igh=n.
    integer, intent(in) :: igh  !! low and igh are integers determined by the balancing
                                !! subroutine  balanc.  if  balanc  has not been used,
                                !! set low=1, igh=n.
    real(wp), intent(inout) :: h(nm, n) !! On input: contains the upper hessenberg matrix.  information about
                                        !! the transformations used in the reduction to hessenberg
                                        !! form by  elmhes  or  orthes, if performed, is stored
                                        !! in the remaining triangle under the hessenberg matrix.
                                        !!
                                        !! On output: has been destroyed.  therefore, it must be saved
                                        !! before calling `hqr` if subsequent calculation and
                                        !! back transformation of eigenvectors is to be performed.
    real(wp), intent(out) :: wr(n) !! the real parts of the eigenvalues.  the eigenvalues
                                   !! are unordered except that complex conjugate pairs
                                   !! of values appear consecutively with the eigenvalue
                                   !! having the positive imaginary part first.  if an
                                   !! error exit is made, the eigenvalues should be correct
                                   !! for indices ierr+1,...,n.
    real(wp), intent(out) :: wi(n) !! the imaginary parts of the eigenvalues.  the eigenvalues
                                   !! are unordered except that complex conjugate pairs
                                   !! of values appear consecutively with the eigenvalue
                                   !! having the positive imaginary part first.  if an
                                   !! error exit is made, the eigenvalues should be correct
                                   !! for indices ierr+1,...,n.
    integer, intent(out) :: ierr !! is set to:
                                 !!
                                 !!  * zero -- for normal return,
                                 !!  * j -- if the limit of 30*n iterations is exhausted
                                 !!    while the j-th eigenvalue is being sought.

    integer :: i, j, k, l, m, en, ll, mm, na, &
               itn, its, mp2, enm2
    real(wp) :: p, q, r, s, t, w, x, y, zz, norm, &
                tst1, tst2
    logical :: notlas

    ierr = 0
    norm = 0.0_wp
    k = 1

    ! store roots isolated by balance and compute matrix norm
    do i = 1, n
        do j = k, n
            norm = norm + abs(h(i, j))
        end do
        k = i
        if (i < low .or. i > igh) then
            wr(i) = h(i, i)
            wi(i) = 0.0_wp
        end if
    end do

    en = igh
    t = 0.0_wp
    itn = 30*n

    do
        ! search for next eigenvalues
        if (en < low) return
        its = 0
        na = en - 1
        enm2 = na - 1
        do
            ! look for single small sub-diagonal element
            ! for l=en step -1 until low do --
            do ll = low, en
                l = en + low - ll
                if (l == low) exit
                s = abs(h(l - 1, l - 1)) + abs(h(l, l))
                if (s == 0.0_wp) s = norm
                tst1 = s
                tst2 = tst1 + abs(h(l, l - 1))
                if (tst2 == tst1) exit
            end do
            ! form shift
            x = h(en, en)
            if (l == en) then
                ! one root found
                wr(en) = x + t
                wi(en) = 0.0_wp
                en = na
            else
                y = h(na, na)
                w = h(en, na)*h(na, en)
                if (l == na) then
                    ! two roots found
                    p = (y - x)/2.0_wp
                    q = p*p + w
                    zz = sqrt(abs(q))
                    x = x + t
                    if (q < 0.0_wp) then
                        ! complex pair
                        wr(na) = x + p
                        wr(en) = x + p
                        wi(na) = zz
                        wi(en) = -zz
                    else
                        ! real pair
                        zz = p + sign(zz, p)
                        wr(na) = x + zz
                        wr(en) = wr(na)
                        if (zz /= 0.0_wp) wr(en) = x - w/zz
                        wi(na) = 0.0_wp
                        wi(en) = 0.0_wp
                    end if
                    en = enm2
                elseif (itn == 0) then
                    ! set error -- all eigenvalues have not
                    ! converged after 30*n iterations
                    ierr = en
                    return
                else
                    if (its == 10 .or. its == 20) then
                        ! form exceptional shift
                        t = t + x
                        do i = low, en
                            h(i, i) = h(i, i) - x
                        end do
                        s = abs(h(en, na)) + abs(h(na, enm2))
                        x = 0.75_wp*s
                        y = x
                        w = -0.4375_wp*s*s
                    end if
                    its = its + 1
                    itn = itn - 1
                    ! look for two consecutive small
                    ! sub-diagonal elements.
                    ! for m=en-2 step -1 until l do --
                    do mm = l, enm2
                        m = enm2 + l - mm
                        zz = h(m, m)
                        r = x - zz
                        s = y - zz
                        p = (r*s - w)/h(m + 1, m) + h(m, m + 1)
                        q = h(m + 1, m + 1) - zz - r - s
                        r = h(m + 2, m + 1)
                        s = abs(p) + abs(q) + abs(r)
                        p = p/s
                        q = q/s
                        r = r/s
                        if (m == l) exit
                        tst1 = abs(p)*(abs(h(m - 1, m - 1)) + abs(zz) + abs(h(m + 1, m + 1)))
                        tst2 = tst1 + abs(h(m, m - 1))*(abs(q) + abs(r))
                        if (tst2 == tst1) exit
                    end do

                    mp2 = m + 2

                    do i = mp2, en
                        h(i, i - 2) = 0.0_wp
                        if (i /= mp2) h(i, i - 3) = 0.0_wp
                    end do
                    ! double qr step involving rows l to en and
                    ! columns m to en
                    do k = m, na
                        notlas = k /= na
                        if (k /= m) then
                            p = h(k, k - 1)
                            q = h(k + 1, k - 1)
                            r = 0.0_wp
                            if (notlas) r = h(k + 2, k - 1)
                            x = abs(p) + abs(q) + abs(r)
                            if (x == 0.0_wp) cycle
                            p = p/x
                            q = q/x
                            r = r/x
                        end if
                        s = sign(sqrt(p*p + q*q + r*r), p)
                        if (k == m) then
                            if (l /= m) h(k, k - 1) = -h(k, k - 1)
                        else
                            h(k, k - 1) = -s*x
                        end if
                        p = p + s
                        x = p/s
                        y = q/s
                        zz = r/s
                        q = q/p
                        r = r/p
                        if (notlas) then
                            ! row modification
                            do j = k, en
                                p = h(k, j) + q*h(k + 1, j) + r*h(k + 2, j)
                                h(k, j) = h(k, j) - p*x
                                h(k + 1, j) = h(k + 1, j) - p*y
                                h(k + 2, j) = h(k + 2, j) - p*zz
                            end do
                            j = min(en, k + 3)
                            ! column modification
                            do i = l, j
                                p = x*h(i, k) + y*h(i, k + 1) + zz*h(i, k + 2)
                                h(i, k) = h(i, k) - p
                                h(i, k + 1) = h(i, k + 1) - p*q
                                h(i, k + 2) = h(i, k + 2) - p*r
                            end do
                        else
                            ! row modification
                            do j = k, en
                                p = h(k, j) + q*h(k + 1, j)
                                h(k, j) = h(k, j) - p*x
                                h(k + 1, j) = h(k + 1, j) - p*y
                            end do
                            j = min(en, k + 3)
                            ! column modification
                            do i = l, j
                                p = x*h(i, k) + y*h(i, k + 1)
                                h(i, k) = h(i, k) - p
                                h(i, k + 1) = h(i, k + 1) - p*q
                            end do
                        end if
                    end do
                    cycle
                end if
            end if
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

    integer, intent(in) :: n !! polynomial degree
    real(wp), dimension(n+1), intent(in) :: p !! polynomial coefficients array, in order of decreasing powers
    real(wp), dimension(n), intent(out) :: wr !! real parts of roots
    real(wp), dimension(n), intent(out) :: wi !! imaginary parts of roots
    integer, intent(out) :: info !! output from the lapack solver.
                                 !!
                                 !! * if `info=0` the routine converged.
                                 !! * if `info=-999`, then the leading coefficient is zero.

    integer :: i !! counter
    real(wp), allocatable, dimension(:,:) :: a !! companion matrix
    real(wp), allocatable, dimension(:) :: work !! work array
    real(wp), dimension(1) :: vl, vr !! not used here

#ifdef REAL32
    interface
        subroutine sgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
            implicit none
            character :: jobvl, jobvr
            integer :: info, lda, ldvl, ldvr, lwork, n
            real :: a(lda, *), vl(ldvl, *), vr(ldvr, *), wi(*), work(*), wr(*)
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
            double precision :: a(lda, *), vl(ldvl, *), vr(ldvr, *), wi(*), work(*), wr(*)
        end subroutine dgeev
    end interface
#endif

    ! error check:
    if (abs(p(1)) == 0.0_wp) then
        info = -999
        return
    end if

    ! allocate the work array:
    allocate (work(3*n))

    ! create the companion matrix
    allocate (a(n, n))
    a = 0.0_wp
    do i = 1, n - 1
        a(i, i + 1) = 1.0_wp
    end do
    do i = n, 1, -1
        a(n, n - i + 1) = -p(i + 1)/p(1)
    end do

    ! call the lapack solver:
#ifdef REAL32
    ! single precision
    call sgeev('N', 'N', n, a, n, wr, wi, vl, 1, vr, 1, work, 3*n, info)
#elif REAL128
    error stop 'do not have a quad solver in lapack'
#else
    ! by default, use double precision:
    call dgeev('N', 'N', n, a, n, wr, wi, vl, 1, vr, 1, work, 3*n, info)
#endif

end subroutine polyroots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve for the roots of a complex polynomial equation by
!  computing the eigenvalues of the companion matrix.
!
!  This one uses LAPACK for the eigen solver (`cgeev` or `zgeev`).
!
!### Reference
!  * Based on [[polyroots]]
!
!### History
!  * Jacob Williams, 9/22/2022 : created this routine.
!
!@note Works only for single and double precision.

subroutine cpolyroots(n, p, w, info)

    implicit none

    integer, intent(in) :: n !! polynomial degree
    complex(wp), dimension(n+1), intent(in) :: p !! polynomial coefficients array, in order of decreasing powers
    complex(wp), dimension(n), intent(out) :: w  !! roots
    integer, intent(out) :: info !! output from the lapack solver.
                                 !!
                                 !! * if `info=0` the routine converged.
                                 !! * if `info=-999`, then the leading coefficient is zero.

    integer :: i !! counter
    complex(wp), allocatable, dimension(:,:) :: a !! companion matrix
    complex(wp), allocatable, dimension(:) :: work !! work array
    real(wp), allocatable, dimension(:) :: rwork !! rwork array (2*n)
    complex(wp), dimension(1) :: vl, vr !! not used here

#ifdef REAL32
    interface
        subroutine cgeev( jobvl, jobvr, n, a, lda, w,      vl, ldvl, vr, ldvr, work, lwork, rwork, info )
            implicit none
            character :: jobvl, jobvr
            integer :: info, lda, ldvl, ldvr, lwork, n
            real :: rwork( * )
            complex :: a( lda, * ), vl( ldvl, * ), vr( ldvr, * ), w( * ), work( * )
        end subroutine cgeev
    end interface
#elif REAL128
    ! do not have a quad solver in lapack
#else
    interface
        subroutine zgeev( jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info )
            implicit none
            character :: jobvl, jobvr
            integer :: info, lda, ldvl, ldvr, lwork, n
            double precision :: rwork( * )
            complex*16 :: a( lda, * ), vl( ldvl, * ), vr( ldvr, * ), w( * ), work( * )
        end subroutine zgeev
    end interface
#endif

    ! error check:
    if (abs(p(1)) == 0.0_wp) then
        info = -999
        return
    end if

    ! allocate the work array:
    allocate (work(3*n))
    allocate(rwork(2*n))

    ! create the companion matrix
    allocate (a(n, n))
    a = 0.0_wp
    do i = 1, n - 1
        a(i, i + 1) = 1.0_wp
    end do
    do i = n, 1, -1
        a(n, n - i + 1) = -p(i + 1)/p(1)
    end do

    ! call the lapack solver:
#ifdef REAL32
    ! single precision
    call cgeev('N', 'N', n, a, n, w, vl, 1, vr, 1, work, 3*n, rwork, info)
#elif REAL128
    error stop 'do not have a quad solver in lapack'
#else
    ! by default, use double precision:
    call zgeev('N', 'N', n, a, n, w, vl, 1, vr, 1, work, 3*n, rwork, info)
#endif

end subroutine cpolyroots
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

subroutine cpqr79(ndeg, coeff, root, ierr)
    implicit none

    integer, intent(in) :: ndeg !! degree of polynomial
    complex(wp), dimension(ndeg+1), intent(in) :: coeff !! `(NDEG+1)` coefficients in descending order.  i.e.,
                                                        !! `P(Z)= COEFF(1)*(Z**NDEG) + COEFF(NDEG)*Z + COEFF(NDEG+1)`
    complex(wp), dimension(ndeg), intent(out) :: root !! `(NDEG)` vector of roots
    integer, intent(out) :: ierr !! Output Error Code.
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

    complex(wp) :: scale, c
    integer :: k, khr, khi, kwr, kwi, kad, kj, km1
    real(wp), dimension(:), allocatable :: work !! work array of dimension `2*NDEG*(NDEG+1)`

    ierr = 0
    if (abs(coeff(1)) == 0.0_wp) then
        ierr = 2
        write (*, *) 'leading coefficient is zero.'
        return
    end if

    if (ndeg <= 0) then
        ierr = 3
        write (*, *) 'degree invalid.'
        return
    end if

    if (ndeg == 1) then
        root(1) = -coeff(2)/coeff(1)
        return
    end if

    ! allocate work array:
    allocate (work(2*NDEG*(NDEG + 1)))

    scale = 1.0_wp/coeff(1)
    khr = 1
    khi = khr + ndeg*ndeg
    kwr = khi + khi - khr
    kwi = kwr + ndeg

    do k = 1, kwr
        work(k) = 0.0_wp
    end do

    do k = 1, ndeg
        kad = (k - 1)*ndeg + 1
        c = scale*coeff(k + 1)
        work(kad) = -real(c, wp)
        kj = khi + kad - 1
        work(kj) = -aimag(c)
        if (k /= ndeg) work(kad + k) = 1.0_wp
    end do

    call comqr(ndeg, ndeg, 1, ndeg, work(khr), work(khi), work(kwr), work(kwi), ierr)

    if (ierr /= 0) then
        write (*, *) 'no convergence in 30 qr iterations. ierr = ', ierr
        ierr = 1
        return
    end if

    do k = 1, ndeg
        km1 = k - 1
        root(k) = cmplx(work(kwr + km1), work(kwi + km1), wp)
    end do

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
!  * calls [[pythag]] for sqrt(a*a + b*b) .
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
!  * From EISPACK. this version dated august 1983.
!    questions and comments should be directed to burton s. garbow,
!    mathematics and computer science div, argonne national laboratory
!  * Jacob Williams, 9/14/2022 : modernized this code

subroutine comqr(nm, n, low, igh, hr, hi, wr, wi, ierr)
    implicit none

    integer, intent(in) :: nm !! row dimension of two-dimensional array parameters
    integer, intent(in) :: n !! the order of the matrix
    integer, intent(in) :: low !! integer determined by the balancing
                               !! subroutine  cbal.  if  cbal  has not been used,
                               !! set low=1
    integer, intent(in) :: igh !! integer determined by the balancing
                               !! subroutine  cbal.  if  cbal  has not been used,
                               !! igh=n.
    real(wp), intent(inout) :: hr(nm, n) !! On input: hr and hi contain the real and imaginary parts,
                                         !! respectively, of the complex upper hessenberg matrix.
                                         !! their lower triangles below the subdiagonal contain
                                         !! information about the unitary transformations used in
                                         !! the reduction by  corth, if performed.
                                         !!
                                         !! On Output: the upper hessenberg portions of hr and hi have been
                                         !! destroyed.  therefore, they must be saved before
                                         !! calling  comqr  if subsequent calculation of
                                         !! eigenvectors is to be performed.
    real(wp), intent(inout) :: hi(nm, n) !! See `hr` description
    real(wp), intent(out) :: wr(n) !! the real parts of the eigenvalues.  if an error
                                   !! exit is made, the eigenvalues should be correct
                                   !! for indices `ierr+1,...,n`.
    real(wp), intent(out) :: wi(n) !! the imaginary parts of the eigenvalues.  if an error
                                   !! exit is made, the eigenvalues should be correct
                                   !! for indices `ierr+1,...,n`.
    integer, intent(out) :: ierr !! is set to:
                                 !!
                                 !!  * 0 -- for normal return
                                 !!  * j -- if the limit of 30*n iterations is exhausted
                                 !!    while the j-th eigenvalue is being sought.

    integer :: i, j, l, en, ll, itn, its, lp1, enm1
    real(wp) :: si, sr, ti, tr, xi, xr, yi, yr, zzi, &
                zzr, norm, tst1, tst2, xr2, xi2

    ierr = 0
    if (low /= igh) then
        ! create real subdiagonal elements
        l = low + 1
        do i = l, igh
            ll = min(i + 1, igh)
            if (hi(i, i - 1) /= 0.0_wp) then
                norm = pythag(hr(i, i - 1), hi(i, i - 1))
                yr = hr(i, i - 1)/norm
                yi = hi(i, i - 1)/norm
                hr(i, i - 1) = norm
                hi(i, i - 1) = 0.0_wp
                do j = i, igh
                    si = yr*hi(i, j) - yi*hr(i, j)
                    hr(i, j) = yr*hr(i, j) + yi*hi(i, j)
                    hi(i, j) = si
                end do
                do j = low, ll
                    si = yr*hi(j, i) + yi*hr(j, i)
                    hr(j, i) = yr*hr(j, i) - yi*hi(j, i)
                    hi(j, i) = si
                end do
            end if
        end do
    end if
    ! store roots isolated by cbal
    do i = 1, n
        if (i < low .or. i > igh) then
            wr(i) = hr(i, i)
            wi(i) = hi(i, i)
        end if
    end do

    en = igh
    tr = 0.0_wp
    ti = 0.0_wp
    itn = 30*n

    main: do

        if (en < low) return

        ! search for next eigenvalue
        its = 0
        enm1 = en - 1

        do

            ! look for single small sub-diagonal element
            ! for l=en step -1 until low d0 --
            do ll = low, en
                l = en + low - ll
                if (l == low) exit
                tst1 = abs(hr(l - 1, l - 1)) + abs(hi(l - 1, l - 1)) + abs(hr(l, l)) + abs(hi(l, l))
                tst2 = tst1 + abs(hr(l, l - 1))
                if (tst2 == tst1) exit
            end do

            ! form shift
            if (l == en) then
                ! a root found
                wr(en) = hr(en, en) + tr
                wi(en) = hi(en, en) + ti
                en = enm1
                cycle main
            elseif (itn == 0) then
                ! set error -- all eigenvalues have not converged after 30*n iterations
                ierr = en
                return
            else
                if (its == 10 .or. its == 20) then
                    ! form exceptional shift
                    sr = abs(hr(en, enm1)) + abs(hr(enm1, en - 2))
                    si = 0.0_wp
                else
                    sr = hr(en, en)
                    si = hi(en, en)
                    xr = hr(enm1, en)*hr(en, enm1)
                    xi = hi(enm1, en)*hr(en, enm1)
                    if (xr /= 0.0_wp .or. xi /= 0.0_wp) then
                        yr = (hr(enm1, enm1) - sr)/2.0_wp
                        yi = (hi(enm1, enm1) - si)/2.0_wp
                        call csroot(yr**2 - yi**2 + xr, 2.0_wp*yr*yi + xi, zzr, zzi)
                        if (yr*zzr + yi*zzi < 0.0_wp) then
                            zzr = -zzr
                            zzi = -zzi
                        end if
                        call cdiv(xr, xi, yr + zzr, yi + zzi, xr2, xi2)
                        xr = xr2
                        xi = xi2
                        sr = sr - xr
                        si = si - xi
                    end if
                end if

                do i = low, en
                    hr(i, i) = hr(i, i) - sr
                    hi(i, i) = hi(i, i) - si
                end do

                tr = tr + sr
                ti = ti + si
                its = its + 1
                itn = itn - 1
                ! reduce to triangle (rows)
                lp1 = l + 1

                do i = lp1, en
                    sr = hr(i, i - 1)
                    hr(i, i - 1) = 0.0_wp
                    norm = pythag(pythag(hr(i - 1, i - 1), hi(i - 1, i - 1)), sr)
                    xr = hr(i - 1, i - 1)/norm
                    wr(i - 1) = xr
                    xi = hi(i - 1, i - 1)/norm
                    wi(i - 1) = xi
                    hr(i - 1, i - 1) = norm
                    hi(i - 1, i - 1) = 0.0_wp
                    hi(i, i - 1) = sr/norm

                    do j = i, en
                        yr = hr(i - 1, j)
                        yi = hi(i - 1, j)
                        zzr = hr(i, j)
                        zzi = hi(i, j)
                        hr(i - 1, j) = xr*yr + xi*yi + hi(i, i - 1)*zzr
                        hi(i - 1, j) = xr*yi - xi*yr + hi(i, i - 1)*zzi
                        hr(i, j) = xr*zzr - xi*zzi - hi(i, i - 1)*yr
                        hi(i, j) = xr*zzi + xi*zzr - hi(i, i - 1)*yi
                    end do

                end do

                si = hi(en, en)
                if (si /= 0.0_wp) then
                    norm = pythag(hr(en, en), si)
                    sr = hr(en, en)/norm
                    si = si/norm
                    hr(en, en) = norm
                    hi(en, en) = 0.0_wp
                end if
                ! inverse operation (columns)
                do j = lp1, en
                    xr = wr(j - 1)
                    xi = wi(j - 1)

                    do i = l, j
                        yr = hr(i, j - 1)
                        yi = 0.0_wp
                        zzr = hr(i, j)
                        zzi = hi(i, j)
                        if (i /= j) then
                            yi = hi(i, j - 1)
                            hi(i, j - 1) = xr*yi + xi*yr + hi(j, j - 1)*zzi
                        end if
                        hr(i, j - 1) = xr*yr - xi*yi + hi(j, j - 1)*zzr
                        hr(i, j) = xr*zzr + xi*zzi - hi(j, j - 1)*yr
                        hi(i, j) = xr*zzi - xi*zzr - hi(j, j - 1)*yi
                    end do

                end do

                if (si /= 0.0_wp) then
                    do i = l, en
                        yr = hr(i, en)
                        yi = hi(i, en)
                        hr(i, en) = sr*yr - si*yi
                        hi(i, en) = sr*yi + si*yr
                    end do
                end if
            end if

        end do

    end do main

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

pure subroutine csroot(xr, xi, yr, yi)
    implicit none

    real(wp), intent(in) :: xr, xi
    real(wp), intent(out) :: yr, yi

    real(wp) :: s, tr, ti

    ! branch chosen so that yr >= 0.0 and sign(yi) == sign(xi)
    tr = xr
    ti = xi
    s = sqrt(0.5_wp*(pythag(tr, ti) + abs(tr)))
    if (tr >= 0.0_wp) yr = s
    if (ti < 0.0_wp) s = -s
    if (tr <= 0.0_wp) yi = s
    if (tr < 0.0_wp) yr = 0.5_wp*(ti/yi)
    if (tr > 0.0_wp) yi = 0.5_wp*(ti/yr)

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

pure subroutine cdiv(ar, ai, br, bi, cr, ci)
    implicit none

    real(wp), intent(in) :: ar, ai, br, bi
    real(wp), intent(out) :: cr, ci

    real(wp) :: s, ars, ais, brs, bis

    s = abs(br) + abs(bi)
    ars = ar/s
    ais = ai/s
    brs = br/s
    bis = bi/s
    s = brs**2 + bis**2
    cr = (ars*brs + ais*bis)/s
    ci = (ais*brs - ars*bis)/s

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

pure real(wp) function pythag(a, b)
    implicit none

    real(wp), intent(in) :: a, b

    real(wp) :: p, q, r, s, t

    p = max(abs(a), abs(b))
    q = min(abs(a), abs(b))

    if (q /= 0.0_wp) then
        do
            r = (q/p)**2
            t = 4.0_wp + r
            if (t == 4.0_wp) exit
            s = r/t
            p = p + 2.0_wp*p*s
            q = q*s
        end do
    end if
    pythag = p

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

subroutine newton_root_polish_real(n, p, zr, zi, ftol, ztol, maxiter, istat)

    implicit none

    integer, intent(in) :: n                  !! degree of polynomial
    real(wp), dimension(n+1), intent(in) :: p !! vector of coefficients in order of decreasing powers
    real(wp), intent(inout) :: zr             !! real part of the zero
    real(wp), intent(inout) :: zi             !! imaginary part of the zero
    real(wp), intent(in) :: ftol              !! convergence tolerance for the root
    real(wp), intent(in) :: ztol              !! convergence tolerance for `x`
    integer, intent(in) :: maxiter            !! maximum number of iterations
    integer, intent(out) :: istat             !! status flag:
                                              !!
                                              !! * 0  = converged in `f`
                                              !! * 1  = converged in `x`
                                              !! * -1 = singular
                                              !! * -2 = max iterations reached

    complex(wp) :: z !! complex number for `(zr,zi)`
    complex(wp) :: f !! function evaluation
    complex(wp) :: z_prev !! previous value of `z`
    complex(wp) :: z_best !! best `z` so far
    complex(wp) :: f_best !! best `f` so far
    complex(wp) :: dfdx !! derivative evaluation
    integer :: i !! counter

    real(wp), parameter :: alpha = 1.0_wp !! newton step size

    ! first evaluate initial point:
    z = cmplx(zr, zi, wp)
    call func(z, f, dfdx)

    ! initialize:
    istat = 0
    z_prev = z
    z_best = z
    f_best = f

    main: do i = 1, maxiter

        if (i > 1) call func(z, f, dfdx)
        if (abs(f) < abs(f_best)) then
            ! best so far
            zr = real(z, wp)
            zi = aimag(z)
            z_best = z
            f_best = f
        end if
        if (abs(f) <= ftol) exit main

        if (i == maxiter) then ! max iterations reached
            istat = -2
            exit main
        end if

        if (dfdx == 0.0_wp) then  ! can't proceed
            istat = -1
            exit main
        end if

        ! Newton correction for next step:
        z = z - alpha*(f/dfdx)

        if (abs(z - z_prev) <= ztol) then
            ! convergence in x. try this point and see if there is any improvement
            istat = 1
            call func(z, f, dfdx)
            if (abs(f) < abs(f_best)) then ! best so far
                zr = real(z, wp)
                zi = aimag(z)
            end if
            exit main
        end if
        z_prev = z

    end do main

contains

    subroutine func(x, f, dfdx)

        !! evaluate the polynomial:
        !! `f = p(1)*x**n + p(2)*x**(n-1) + ... + p(n)*x + p(n+1)`
        !! and its derivative using Horner's Rule.
        !!
        !! See: "Roundoff in polynomial evaluation", W. Kahan, 1986
        !! https://rosettacode.org/wiki/Horner%27s_rule_for_polynomial_evaluation#Fortran

        implicit none

        complex(wp), intent(in) :: x
        complex(wp), intent(out) :: f    !! function value at `x`
        complex(wp), intent(out) :: dfdx !! function derivative at `x`

        integer :: i !! counter

        f = p(1)
        dfdx = 0.0_wp
        do i = 2, n + 1
            dfdx = dfdx*x + f
            f = f*x + p(i)
        end do

    end subroutine func

end subroutine newton_root_polish_real
!*****************************************************************************************

!*****************************************************************************************
!>
!  "Polish" a root using a complex Newton Raphson method.
!  This routine will perform a Newton iteration and update the roots only if they improve,
!  otherwise, they are left as is.
!
!@note Same as [[newton_root_polish_real]] except that `p` is `complex(wp)`

subroutine newton_root_polish_complex(n, p, zr, zi, ftol, ztol, maxiter, istat)

    implicit none

    integer, intent(in) :: n                     !! degree of polynomial
    complex(wp), dimension(n+1), intent(in) :: p !! vector of coefficients in order of decreasing powers
    real(wp), intent(inout) :: zr                !! real part of the zero
    real(wp), intent(inout) :: zi                !! imaginary part of the zero
    real(wp), intent(in) :: ftol                 !! convergence tolerance for the root
    real(wp), intent(in) :: ztol                 !! convergence tolerance for `x`
    integer, intent(in) :: maxiter               !! maximum number of iterations
    integer, intent(out) :: istat                !! status flag:
                                                 !!
                                                 !! * 0  = converged in `f`
                                                 !! * 1  = converged in `x`
                                                 !! * -1 = singular
                                                 !! * -2 = max iterations reached

    complex(wp) :: z !! complex number for `(zr,zi)`
    complex(wp) :: f !! function evaluation
    complex(wp) :: z_prev !! previous value of `z`
    complex(wp) :: z_best !! best `z` so far
    complex(wp) :: f_best !! best `f` so far
    complex(wp) :: dfdx !! derivative evaluation
    integer :: i !! counter

    real(wp), parameter :: alpha = 1.0_wp !! newton step size

    ! first evaluate initial point:
    z = cmplx(zr, zi, wp)
    call func(z, f, dfdx)

    ! initialize:
    istat = 0
    z_prev = z
    z_best = z
    f_best = f

    main: do i = 1, maxiter

        if (i > 1) call func(z, f, dfdx)
        if (abs(f) < abs(f_best)) then
            ! best so far
            zr = real(z, wp)
            zi = aimag(z)
            z_best = z
            f_best = f
        end if
        if (abs(f) <= ftol) exit main

        if (i == maxiter) then ! max iterations reached
            istat = -2
            exit main
        end if

        if (dfdx == 0.0_wp) then  ! can't proceed
            istat = -1
            exit main
        end if

        ! Newton correction for next step:
        z = z - alpha*(f/dfdx)

        if (abs(z - z_prev) <= ztol) then
            ! convergence in x. try this point and see if there is any improvement
            istat = 1
            call func(z, f, dfdx)
            if (abs(f) < abs(f_best)) then ! best so far
                zr = real(z, wp)
                zi = aimag(z)
            end if
            exit main
        end if
        z_prev = z

    end do main

contains

    subroutine func(x, f, dfdx)

        !! evaluate the polynomial:
        !! `f = p(1)*x**n + p(2)*x**(n-1) + ... + p(n)*x + p(n+1)`
        !! and its derivative using Horner's Rule.
        !!
        !! See: "Roundoff in polynomial evaluation", W. Kahan, 1986
        !! https://rosettacode.org/wiki/Horner%27s_rule_for_polynomial_evaluation#Fortran

        implicit none

        complex(wp), intent(in) :: x
        complex(wp), intent(out) :: f    !! function value at `x`
        complex(wp), intent(out) :: dfdx !! function derivative at `x`

        integer :: i !! counter

        f = p(1)
        dfdx = 0.0_wp
        do i = 2, n + 1
            dfdx = dfdx*x + f
            f = f*x + p(i)
        end do

    end subroutine func

end subroutine newton_root_polish_complex
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subroutine finds roots of a complex polynomial.
!  It uses a new dynamic root finding algorithm (see the Paper).
!
!  It can use Laguerre's method (subroutine [[cmplx_laguerre]])
!  or Laguerre->SG->Newton method (subroutine
!  [[cmplx_laguerre2newton]] - this is default choice) to find
!  roots. It divides polynomial one by one by found roots. At the
!  end it finds last root from Viete's formula for quadratic
!  equation. Finally, it polishes all found roots using a full
!  polynomial and Newton's or Laguerre's method (default is
!  Laguerre's - subroutine [[cmplx_laguerre]]).
!  You can change default choices by commenting out and uncommenting
!  certain lines in the code below.
!
!### Reference
!  * J. Skowron & A. Gould,
!    "[General Complex Polynomial Root Solver and Its Further Optimization for Binary Microlenses](https://arxiv.org/pdf/1203.1034.pdf)"
!    (2012)
!
!### History
!  * Original code here (Apache license): http://www.astrouw.edu.pl/~jskowron/cmplx_roots_sg/
!  * Jacob Williams, 9/18/2022 : refactored this code a bit
!
!### Notes:
!
! * we solve for the last root with Viete's formula rather
!   than doing full Laguerre step (which is time consuming
!   and unnecessary)
! * we do not introduce any preference to real roots
! * in Laguerre implementation we omit unneccesarry calculation of
!   absolute values of denominator
! * we do not sort roots.

subroutine cmplx_roots_gen(degree, poly, roots, polish_roots_after, use_roots_as_starting_points)

    implicit none

    integer, intent(in) :: degree !! degree of the polynomial and size of 'roots' array
    complex(wp), dimension(degree+1), intent(in) :: poly !! coeffs of the polynomial, in order of increasing powers.
    complex(wp), dimension(degree), intent(inout) :: roots !! array which will hold all roots that had been found.
                                                           !! If the flag 'use_roots_as_starting_points' is set to
                                                           !! .true., then instead of point (0,0) we use value from
                                                           !! this array as starting point for cmplx_laguerre
    logical, intent(in), optional :: polish_roots_after !! after all roots have been found by dividing
                                                        !! original polynomial by each root found,
                                                        !! you can opt in to polish all roots using full
                                                        !! polynomial. [default is false]
    logical, intent(in), optional :: use_roots_as_starting_points !! usually we start Laguerre's
                                                                  !! method from point (0,0), but you can decide to use the
                                                                  !! values of 'roots' array as starting point for each new
                                                                  !! root that is searched for. This is useful if you have
                                                                  !! very rough idea where some of the roots can be.
                                                                  !! [default is false]

    complex(wp), dimension(:), allocatable :: poly2 !! `degree+1` array
    integer :: i, n, iter
    logical :: success
    complex(wp) :: coef, prev

    integer, parameter :: MAX_ITERS=50
    ! constants needed to break cycles in the scheme
    integer, parameter :: FRAC_JUMP_EVERY=10
    integer, parameter :: FRAC_JUMP_LEN=10
    real(wp), dimension(FRAC_JUMP_LEN), parameter :: FRAC_JUMPS=&
        [0.64109297_wp, 0.91577881_wp, 0.25921289_wp,  0.50487203_wp, 0.08177045_wp, &
         0.13653241_wp,  0.306162_wp , 0.37794326_wp, 0.04618805_wp,  0.75132137_wp] !! some random numbers
    real(wp), parameter :: FRAC_ERR = 10.0_wp * eps  !! fractional error
                                                     !! (see. Adams 1967 Eqs 9 and 10)
                                                     !! [2.0d-15 in original code]
    complex(wp), parameter :: zero = cmplx(0.0_wp,0.0_wp,wp)
    complex(wp), parameter :: c_one=cmplx(1.0_wp,0.0_wp,wp)

    ! initialize starting points
    if (present(use_roots_as_starting_points)) then
        if (.not.use_roots_as_starting_points) roots = zero
    else
        roots = zero
    end if

    ! skip small degree polynomials from doing Laguerre's method
    if (degree<=1) then
      if (degree==1) roots(1)=-poly(1)/poly(2)
      return
    endif

    allocate(poly2(degree+1))
    poly2=poly

    do n=degree, 3, -1

      ! find root with Laguerre's method
      !call cmplx_laguerre(poly2, n, roots(n), iter, success)
      ! or
      ! find root with (Laguerre's method -> SG method -> Newton's method)
      call cmplx_laguerre2newton(poly2, n, roots(n), iter, success, 2)
      if (.not.success) then
        roots(n)=zero
        call cmplx_laguerre(poly2, n, roots(n), iter, success)
      endif

      ! divide the polynomial by this root
      coef=poly2(n+1)
      do i=n,1,-1
        prev=poly2(i)
        poly2(i)=coef
        coef=prev+roots(n)*coef
      enddo
      ! variable coef now holds a remainder - should be close to 0

    enddo

    ! find all but last root with Laguerre's method
    !call cmplx_laguerre(poly2, 2, roots(2), iter, success)
    ! or
    call cmplx_laguerre2newton(poly2, 2, roots(2), iter, success, 2)
    if (.not.success) then
      call solve_quadratic_eq(roots(2),roots(1),poly2)
    else
      ! calculate last root from Viete's formula
      roots(1)=-(roots(2)+poly2(2)/poly2(3))
    endif

    if (present(polish_roots_after)) then
        if (polish_roots_after) then
            do n=1, degree ! polish roots one-by-one with a full polynomial
                call cmplx_laguerre(poly, degree, roots(n), iter, success)
                !call cmplx_newton_spec(poly, degree, roots(n), iter, success)
            enddo
        endif
    end if

    contains

    recursive subroutine cmplx_laguerre(poly, degree, root, iter, success)

    !  Subroutine finds one root of a complex polynomial using
    !  Laguerre's method. In every loop it calculates simplified
    !  Adams' stopping criterion for the value of the polynomial.
    !
    !  For a summary of the method go to:
    !  http://en.wikipedia.org/wiki/Laguerre's_method

    implicit none

    integer, intent(in) :: degree !! a degree of the polynomial
    complex(wp), dimension(degree+1), intent(in)  :: poly !! an array of polynomial cooefs
                                                          !! length = degree+1, poly(1) is constant
                                                          !!```
                                                          !!        1              2             3
                                                          !!   poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
                                                          !!```
    integer, intent(out) :: iter !! number of iterations performed (the number of polynomial
                                 !! evaluations and stopping criterion evaluation)
    complex(wp), intent(inout) :: root !! input: guess for the value of a root
                                       !! output: a root of the polynomial
                                       !!
                                       !! Uses 'root' value as a starting point (!!!!!)
                                       !! Remember to initialize 'root' to some initial guess or to
                                       !! point (0,0) if you have no prior knowledge.

    logical, intent(out) :: success !! is false if routine reaches maximum number of iterations

    real(wp) :: faq !! jump length
    complex(wp) :: p         !! value of polynomial
    complex(wp) :: dp        !! value of 1st derivative
    complex(wp) :: d2p_half  !! value of 2nd derivative
    integer :: i, k
    logical :: good_to_go
    complex(wp) :: denom, denom_sqrt, dx, newroot, fac_netwon, fac_extra, F_half, c_one_nth
    real(wp) :: ek, absroot, abs2p, one_nth, n_1_nth, two_n_div_n_1, stopping_crit2

    iter=0
    success=.true.

    ! next if-endif block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
    !if (.false.) then ! change false-->true if you would like to use caution about having first coefficient == 0
      if (degree<0) then
        write(*,*) 'Error: cmplx_laguerre: degree<0'
        return
      endif
      if (poly(degree+1)==zero) then
        if (degree==0) return
        call cmplx_laguerre(poly, degree-1, root, iter, success)
        return
      endif
      if (degree<=1) then
        if (degree==0) then  ! we know from previous check than poly(1) not equal zero
          success=.false.
          write(*,*) 'Warning: cmplx_laguerre: degree=0 and poly(1)/=0, no roots'
          return
        else
          root=-poly(1)/poly(2)
          return
        endif
      endif
    !endif
    !  end EXTREME failsafe

    good_to_go=.false.
    one_nth=1.0_wp/degree
    n_1_nth=(degree-1.0_wp)*one_nth
    two_n_div_n_1=2.0_wp/n_1_nth
    c_one_nth=cmplx(one_nth,0.0_wp,wp)

    do i=1,MAX_ITERS
      ! prepare stoping criterion
      ek=abs(poly(degree+1))
      absroot=abs(root)
      ! calculate value of polynomial and its first two derivatives
      p  =poly(degree+1)
      dp =zero
      d2p_half=zero
      do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
        d2p_half=dp + d2p_half*root
        dp =p + dp*root
        p  =poly(k)+p*root    ! b_k
        ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
        ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
        ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
        ! Eq 8.
        ek=absroot*ek+abs(p)
      enddo
      iter=iter+1

      abs2p=real(conjg(p)*p)
      if (abs2p==0.0_wp) return
      stopping_crit2=(FRAC_ERR*ek)**2
      if (abs2p<stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967)
        ! do additional iteration if we are less than 10x from stopping criterion
        if (abs2p<0.01d0*stopping_crit2) then
          return ! return immediately, because we are at very good place
        else
          good_to_go=.true. ! do one iteration more
        endif
      else
        good_to_go=.false.  ! reset if we are outside the zone of the root
      endif

      faq=1.0_wp
      denom=zero
      if (dp/=zero) then
        fac_netwon=p/dp
        fac_extra=d2p_half/dp
        F_half=fac_netwon*fac_extra

        denom_sqrt=sqrt(c_one-two_n_div_n_1*F_half)

        !G=dp/p  ! gradient of ln(p)
        !G2=G*G
        !H=G2-2.0_wp*d2p_half/p  ! second derivative of ln(p)
        !denom_sqrt=sqrt( (degree-1)*(degree*H-G2) )

        ! NEXT LINE PROBABLY CAN BE COMMENTED OUT
        if (real(denom_sqrt, wp)>=0.0_wp) then
          ! real part of a square root is positive for probably all compilers. You can
          ! test this on your compiler and if so, you can omit this check
          denom=c_one_nth+n_1_nth*denom_sqrt
        else
          denom=c_one_nth-n_1_nth*denom_sqrt
        endif
      endif
      if (denom==zero) then !test if demoninators are > 0.0 not to divide by zero
        dx=(absroot+1.0_wp)*exp(cmplx(0.0_wp,FRAC_JUMPS(mod(i,FRAC_JUMP_LEN)+1)*2*pi,wp)) ! make some random jump
      else
        dx=fac_netwon/denom
        !dx=degree/denom
      endif

      newroot=root-dx
      if (newroot==root) return ! nothing changes -> return
      if (good_to_go) then       ! this was jump already after stopping criterion was met
        root=newroot
        return
      endif

      if (mod(i,FRAC_JUMP_EVERY)==0) then ! decide whether to do a jump of modified length (to break cycles)
        faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
        newroot=root-faq*dx ! do jump of some semi-random length (0<faq<1)
      endif
      root=newroot
    enddo
    success=.false.
    ! too many iterations here
  end subroutine cmplx_laguerre

  subroutine solve_quadratic_eq(x0,x1,poly)

      ! Quadratic equation solver for complex polynomial (degree=2)

      implicit none

      complex(wp), intent(out) :: x0, x1
      complex(wp), dimension(*), intent(in) :: poly !! coeffs of the polynomial
                                                    !! an array of polynomial cooefs,
                                                    !! length = degree+1, poly(1) is constant
                                                    !!```
                                                    !!        1              2             3
                                                    !!   poly(1) x^0 + poly(2) x^1 + poly(3) x^2
                                                    !!```
      complex(wp) :: a, b, c, b2, delta

      a=poly(3)
      b=poly(2)
      c=poly(1)
      ! quadratic equation: a z^2 + b z + c = 0

      b2=b*b
      delta=sqrt(b2-4.0_wp*(a*c))
      if ( real(conjg(b)*delta, wp)>=0.0_wp ) then  ! scalar product to decide the sign yielding bigger magnitude
        x0=-0.5_wp*(b+delta)
      else
        x0=-0.5_wp*(b-delta)
      endif
      if (x0==cmplx(0.0_wp,0.0_wp,wp)) then
        x1=cmplx(0.0_wp,0.0_wp,wp)
      else ! Viete's formula
        x1=c/x0
        x0=x0/a
      endif

    end subroutine solve_quadratic_eq

  recursive subroutine cmplx_laguerre2newton(poly, degree, root, iter, success, starting_mode)

    !  Subroutine finds one root of a complex polynomial using
    !  Laguerre's method, Second-order General method and Newton's
    !  method - depending on the value of function F, which is a
    !  combination of second derivative, first derivative and
    !  value of polynomial [F=-(p"*p)/(p'p')].
    !
    !  Subroutine has 3 modes of operation. It starts with mode=2
    !  which is the Laguerre's method, and continues until F
    !  becames F<0.50, at which point, it switches to mode=1,
    !  i.e., SG method (see paper). While in the first two
    !  modes, routine calculates stopping criterion once per every
    !  iteration. Switch to the last mode, Newton's method, (mode=0)
    !  happens when becomes F<0.05. In this mode, routine calculates
    !  stopping criterion only once, at the beginning, under an
    !  assumption that we are already very close to the root.
    !  If there are more than 10 iterations in Newton's mode,
    !  it means that in fact we were far from the root, and
    !  routine goes back to Laguerre's method (mode=2).
    !
    !  For a summary of the method see the paper: Skowron & Gould (2012)

    implicit none

    integer, intent(in) :: degree !! a degree of the polynomial
    complex(wp), dimension(degree+1), intent(in)  :: poly !! is an array of polynomial cooefs
                                                               !! length = degree+1, poly(1) is constant
                                                               !!```
                                                               !!        1              2             3
                                                               !!   poly(1) x^0 + poly(2) x^1 + poly(3) x^2 + ...
                                                               !!```
    complex(wp), intent(inout) :: root !! input: guess for the value of a root
                                            !! output: a root of the polynomial
                                            !!
                                            !! Uses 'root' value as a starting point (!!!!!)
                                            !! Remember to initialize 'root' to some initial guess or to
                                            !! point (0,0) if you have no prior knowledge.
    integer, intent(in) :: starting_mode !! this should be by default = 2. However if you
                                         !! choose to start with SG method put 1 instead.
                                         !! Zero will cause the routine to
                                         !! start with Newton for first 10 iterations, and
                                         !! then go back to mode 2.
    integer, intent(out) :: iter !! number of iterations performed (the number of polynomial
                                 !! evaluations and stopping criterion evaluation)
    logical, intent(out) :: success !! is false if routine reaches maximum number of iterations

    real(wp) :: faq ! jump length
    complex(wp) :: p         ! value of polynomial
    complex(wp) :: dp        ! value of 1st derivative
    complex(wp) :: d2p_half  ! value of 2nd derivative
    integer :: i, j, k
    logical :: good_to_go
    complex(wp) :: denom, denom_sqrt, dx, newroot
    real(wp) :: ek, absroot, abs2p, abs2_F_half
    complex(wp) :: fac_netwon, fac_extra, F_half, c_one_nth
    real(wp) :: one_nth, n_1_nth, two_n_div_n_1
    integer :: mode
    real(wp) :: stopping_crit2

    iter=0
    success=.true.
    stopping_crit2 = 0.0_wp  !  value not important, will be initialized anyway on the first loop (because mod(1,10)==1)

    ! next if-endif block is an EXTREME failsafe, not usually needed, and thus turned off in this version.
    !if (.false.)then ! change false-->true if you would like to use caution about having first coefficient == 0
      if (degree<0) then
          write(*,*) 'Error: cmplx_laguerre2newton: degree<0'
          return
      endif
      if (poly(degree+1)==zero) then
          if (degree==0) return
          call cmplx_laguerre2newton(poly, degree-1, root, iter, success, starting_mode)
          return
      endif
      if (degree<=1) then
          if (degree==0) then  ! we know from previous check than poly(1) not equal zero
            success=.false.
            write(*,*) 'Warning: cmplx_laguerre2newton: degree=0 and poly(1)/=0, no roots'
            return
          else
            root=-poly(1)/poly(2)
            return
          endif
      endif
    !endif
    !  end EXTREME failsafe

    j=1
    good_to_go=.false.
    mode=starting_mode  ! mode=2 full laguerre, mode=1 SG, mode=0 newton

    do ! infinite loop, just to be able to come back from newton, if more than 10 iteration there

      !------------------------------------------------------------- mode 2
      if (mode>=2) then  ! LAGUERRE'S METHOD
          one_nth=1.0_wp/degree
          n_1_nth=(degree-1.0_wp)*one_nth
          two_n_div_n_1=2.0_wp/n_1_nth
          c_one_nth=cmplx(one_nth,0.0_wp,wp)

          do i=1,MAX_ITERS  !
          faq=1.0_wp

          ! prepare stoping criterion
          ek=abs(poly(degree+1))
          absroot=abs(root)
          ! calculate value of polynomial and its first two derivatives
          p  =poly(degree+1)
          dp =zero
          d2p_half=zero
          do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
              d2p_half=dp + d2p_half*root
              dp =p + dp*root
              p  =poly(k)+p*root    ! b_k
              ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
              ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
              ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
              ! Eq 8.
              ek=absroot*ek+abs(p)
          enddo
          abs2p=real(conjg(p)*p, wp) !abs(p)
          iter=iter+1
          if (abs2p==0.0_wp) return

          stopping_crit2=(FRAC_ERR*ek)**2
          if (abs2p<stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967)
              ! do additional iteration if we are less than 10x from stopping criterion
              if (abs2p<0.01_wp*stopping_crit2) then ! ten times better than stopping criterion
              return ! return immediately, because we are at very good place
              else
              good_to_go=.true. ! do one iteration more
              endif
          else
              good_to_go=.false. ! reset if we are outside the zone of the root
          endif

          denom=zero
          if (dp/=zero) then
              fac_netwon=p/dp
              fac_extra=d2p_half/dp
              F_half=fac_netwon*fac_extra

              abs2_F_half=real(conjg(F_half)*F_half, wp)
              if (abs2_F_half<=0.0625_wp) then     ! F<0.50, F/2<0.25
              ! go to SG method
              if (abs2_F_half<=0.000625_wp) then ! F<0.05, F/2<0.025
                  mode=0 ! go to Newton's
              else
                  mode=1 ! go to SG
              endif
              endif

              denom_sqrt=sqrt(c_one-two_n_div_n_1*F_half)

              ! NEXT LINE PROBABLY CAN BE COMMENTED OUT
              if (real(denom_sqrt, wp)>=0.0_wp) then
                ! real part of a square root is positive for probably all compilers. You can
                ! test this on your compiler and if so, you can omit this check
                denom=c_one_nth+n_1_nth*denom_sqrt
              else
                denom=c_one_nth-n_1_nth*denom_sqrt
              endif
          endif
          if (denom==zero) then !test if demoninators are > 0.0 not to divide by zero
              dx=(abs(root)+1.0_wp)*exp(cmplx(0.0_wp,FRAC_JUMPS(mod(i,FRAC_JUMP_LEN)+1)*2*pi,wp)) ! make some random jump
          else
              dx=fac_netwon/denom
          endif

          newroot=root-dx
          if (newroot==root) return ! nothing changes -> return
          if (good_to_go) then       ! this was jump already after stopping criterion was met
              root=newroot
              return
          endif

          if (mode/=2) then
              root=newroot
              j=i+1    ! remember iteration index
              exit     ! go to Newton's or SG
          endif

          if (mod(i,FRAC_JUMP_EVERY)==0) then ! decide whether to do a jump of modified length (to break cycles)
              faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
              newroot=root-faq*dx ! do jump of some semi-random length (0<faq<1)
          endif
          root=newroot
          enddo ! do mode 2

          if (i>=MAX_ITERS) then
            success=.false.
            return
          endif

      endif ! if mode 2

      !------------------------------------------------------------- mode 1
      if (mode==1) then  ! SECOND-ORDER GENERAL METHOD (SG)

          do i=j,MAX_ITERS  !
          faq=1.0_wp

          ! calculate value of polynomial and its first two derivatives
          p  =poly(degree+1)
          dp =zero
          d2p_half=zero
          if (mod(i-j,10)==0) then
              ! prepare stoping criterion
              ek=abs(poly(degree+1))
              absroot=abs(root)
              do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                d2p_half=dp + d2p_half*root
                dp =p + dp*root
                p  =poly(k)+p*root    ! b_k
                ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
                ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
                ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
                ! Eq 8.
                ek=absroot*ek+abs(p)
              enddo
              stopping_crit2=(FRAC_ERR*ek)**2
          else
              do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                d2p_half=dp + d2p_half*root
                dp =p + dp*root
                p  =poly(k)+p*root    ! b_k
              enddo
          endif

          abs2p=real(conjg(p)*p, wp) !abs(p)**2
          iter=iter+1
          if (abs2p==0.0_wp) return

          if (abs2p<stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967)
              if (dp==zero) return
                ! do additional iteration if we are less than 10x from stopping criterion
                if (abs2p<0.01_wp*stopping_crit2) then ! ten times better than stopping criterion
                return ! return immediately, because we are at very good place
              else
                good_to_go=.true. ! do one iteration more
              endif
          else
              good_to_go=.false. ! reset if we are outside the zone of the root
          endif

          if (dp==zero) then !test if demoninators are > 0.0 not to divide by zero
              dx=(abs(root)+1.0_wp)*exp(cmplx(0.0_wp,FRAC_JUMPS(mod(i,FRAC_JUMP_LEN)+1)*2*pi,wp)) ! make some random jump
          else
              fac_netwon=p/dp
              fac_extra=d2p_half/dp
              F_half=fac_netwon*fac_extra

              abs2_F_half=real(conjg(F_half)*F_half, wp)
              if (abs2_F_half<=0.000625_wp) then ! F<0.05, F/2<0.025
                mode=0 ! set Newton's, go there after jump
              endif

              dx=fac_netwon*(c_one+F_half)  ! SG
          endif

          newroot=root-dx
          if (newroot==root) return ! nothing changes -> return
          if (good_to_go) then       ! this was jump already after stopping criterion was met
              root=newroot
              return
          endif

          if (mode/=1) then
              root=newroot
              j=i+1    ! remember iteration number
              exit     ! go to Newton's
          endif

          if (mod(i,FRAC_JUMP_EVERY)==0) then ! decide whether to do a jump of modified length (to break cycles)
              faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
              newroot=root-faq*dx ! do jump of some semi-random length (0<faq<1)
          endif
          root=newroot

          enddo ! do mode 1

          if (i>=MAX_ITERS) then
            success=.false.
            return
          endif

      endif ! if mode 1

      !------------------------------------------------------------- mode 0
      if (mode==0) then  ! NEWTON'S METHOD

          do i=j,j+10  ! do only 10 iterations the most, then go back to full Laguerre's
            faq=1.0_wp

            ! calculate value of polynomial and its first two derivatives
            p  =poly(degree+1)
            dp =zero
            if (i==j) then ! calculate stopping crit only once at the begining
                ! prepare stoping criterion
                ek=abs(poly(degree+1))
                absroot=abs(root)
                do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                    dp =p + dp*root
                    p  =poly(k)+p*root    ! b_k
                    ! Adams, Duane A., 1967, "A stopping criterion for polynomial root finding",
                    ! Communications of the ACM, Volume 10 Issue 10, Oct. 1967, p. 655
                    ! ftp://reports.stanford.edu/pub/cstr/reports/cs/tr/67/55/CS-TR-67-55.pdf
                    ! Eq 8.
                    ek=absroot*ek+abs(p)
                enddo
                stopping_crit2=(FRAC_ERR*ek)**2
            else        !
                do k=degree,1,-1 ! Horner Scheme, see for eg.  Numerical Recipes Sec. 5.3 how to evaluate polynomials and derivatives
                    dp =p + dp*root
                    p  =poly(k)+p*root    ! b_k
                enddo
            endif
            abs2p=real(conjg(p)*p, wp) !abs(p)**2
            iter=iter+1
            if (abs2p==0.0_wp) return

            if (abs2p<stopping_crit2) then ! (simplified a little Eq. 10 of Adams 1967)
                if (dp==zero) return
                ! do additional iteration if we are less than 10x from stopping criterion
                if (abs2p<0.01_wp*stopping_crit2) then ! ten times better than stopping criterion
                    return ! return immediately, because we are at very good place
                else
                    good_to_go=.true. ! do one iteration more
                endif
            else
                good_to_go=.false. ! reset if we are outside the zone of the root
            endif

            if (dp==zero) then ! test if demoninators are > 0.0 not to divide by zero
                dx=(abs(root)+1.0_wp)*exp(cmplx(0.0_wp,FRAC_JUMPS(mod(i,FRAC_JUMP_LEN)+1)*2*pi,wp)) ! make some random jump
            else
                dx=p/dp
            endif

            newroot=root-dx
            if (newroot==root) return ! nothing changes -> return
            if (good_to_go) then
                root=newroot
                return
            endif

            ! this loop is done only 10 times. So skip this check
            !if (mod(i,FRAC_JUMP_EVERY)==0) then ! decide whether to do a jump of modified length (to break cycles)
            !  faq=FRAC_JUMPS(mod(i/FRAC_JUMP_EVERY-1,FRAC_JUMP_LEN)+1)
            !  newroot=root-faq*dx ! do jump of some semi-random length (0<faq<1)
            !endif
            root=newroot

          enddo ! do mode 0 10 times

          if (iter>=MAX_ITERS) then
            ! too many iterations here
            success=.false.
            return
          endif
          mode=2 ! go back to Laguerre's. This happens when we were unable to converge in 10 iterations with Newton's

      endif ! if mode 0

    enddo ! end of infinite loop

    success=.false.

  end subroutine cmplx_laguerre2newton

end subroutine cmplx_roots_gen
!*****************************************************************************************

!*****************************************************************************************
!>
!  Finds the zeros of a complex polynomial.
!
!### Reference
!  * Jenkins & Traub,
!    "[Algorithm 419: Zeros of a complex polynomial](https://netlib.org/toms-2014-06-10/419)"
!    Communications of the ACM, Volume 15, Issue 2, Feb. 1972, pp 97-99.
!  * Added changes from remark on algorithm 419 by david h. withers
!    cacm (march 1974) vol 17 no 3 p. 157]
!
!@note the program has been written to reduce the chance of overflow
!      occurring. if it does occur, there is still a possibility that
!      the zerofinder will work provided the overflowed quantity is
!      replaced by a large number.
!
!### History
!  * Jacob Williams, 9/18/2022 : modern Fortran version of this code.

subroutine cpoly(opr,opi,degree,zeror,zeroi,fail)

    implicit none

    integer,intent(in) :: degree !! degree of polynomial
    real(wp), dimension(degree+1), intent(in) :: opr !! vectors of real parts of the coefficients in
                                                     !! order of decreasing powers.
    real(wp), dimension(degree+1), intent(in) :: opi !! vectors of imaginary parts of the coefficients in
                                                     !! order of decreasing powers.
    real(wp), dimension(degree), intent(out) :: zeror !! real parts of the zeros
    real(wp), dimension(degree), intent(out) :: zeroi !! imaginary parts of the zeros
    logical,intent(out) :: fail !! true only if leading coefficient is zero or if cpoly
                                !! has found fewer than `degree` zeros.

    real(wp) :: sr , si , tr , ti , pvr , pvi, xxx , zr , zi , bnd , xx , yy
    real(wp), dimension(degree+1) :: pr , pi , hr, hi , qpr, qpi, qhr , qhi , shr , shi
    logical :: conv
    integer :: cnt1 , cnt2, i , idnn2 , nn

    real(wp), parameter :: base = radix(1.0_wp)
    real(wp), parameter :: eta = eps
    real(wp), parameter :: infin = huge(1.0_wp)
    real(wp), parameter :: smalno = tiny(1.0_wp)
    real(wp), parameter :: are = eta
    real(wp), parameter :: cosr = cos(94.0_wp*deg2rad)  !! -.069756474
    real(wp), parameter :: sinr = sin(86.0_wp*deg2rad)  !! .99756405
    real(wp), parameter :: mre = 2.0_wp*sqrt(2.0_wp)*eta
    real(wp), parameter :: cos45 = cos(45.0_wp*deg2rad) !! .70710678

    if ( opr(1)==0.0_wp .and. opi(1)==0.0_wp ) then
        ! algorithm fails if the leading coefficient is zero.
        fail = .true.
        return
    end if

    xx = cos45
    yy = -xx
    fail = .false.
    nn = degree + 1

    ! remove the zeros at the origin if any.
    do
        if ( opr(nn)/=0.0_wp .or. opi(nn)/=0.0_wp ) then
            exit
        else
            idnn2 = degree - nn + 2
            zeror(idnn2) = 0.0_wp
            zeroi(idnn2) = 0.0_wp
            nn = nn - 1
        endif
    end do

    ! make a copy of the coefficients.
    do i = 1 , nn
        pr(i) = opr(i)
        pi(i) = opi(i)
        shr(i) = cmod(pr(i),pi(i))
    enddo
    ! scale the polynomial.
    bnd = scale(nn,shr,eta,infin,smalno,base)
    if ( bnd/=1.0_wp ) then
        do i = 1 , nn
            pr(i) = bnd*pr(i)
            pi(i) = bnd*pi(i)
        enddo
    endif

    ! start the algorithm for one zero.
    main : do
        if ( nn>2 ) then
            ! calculate bnd, a lower bound on the modulus of the zeros.
            do i = 1 , nn
                shr(i) = cmod(pr(i),pi(i))
            enddo
            bnd = cauchy(nn,shr,shi)
            ! outer loop to control 2 major passes with different sequences
            ! of shifts.
            do cnt1 = 1 , 2
                ! first stage calculation, no shift.
                call noshft(5)
                ! inner loop to select a shift.
                do cnt2 = 1 , 9
                    ! shift is chosen with modulus bnd and amplitude rotated by
                    ! 94 degrees from the previous shift
                    xxx = cosr*xx - sinr*yy
                    yy = sinr*xx + cosr*yy
                    xx = xxx
                    sr = bnd*xx
                    si = bnd*yy
                    ! second stage calculation, fixed shift.
                    call fxshft(10*cnt2,zr,zi,conv)
                    if ( conv ) then
                        ! the second stage jumps directly to the third stage iteration.
                        ! if successful the zero is stored and the polynomial deflated.
                        idnn2 = degree - nn + 2
                        zeror(idnn2) = zr
                        zeroi(idnn2) = zi
                        nn = nn - 1
                        do i = 1 , nn
                            pr(i) = qpr(i)
                            pi(i) = qpi(i)
                        enddo
                        cycle main
                    endif
                    ! if the iteration is unsuccessful another shift is chosen.
                enddo
                ! if 9 shifts fail, the outer loop is repeated with another
                ! sequence of shifts.
            enddo
            ! the zerofinder has failed on two major passes.
            ! return empty handed.
            fail = .true.
            return
        else
            exit
        endif

    end do main

    ! calculate the final zero and return.
    call cdivid(-pr(2),-pi(2),pr(1),pi(1),zeror(degree),zeroi(degree))

contains

subroutine noshft(l1)

    ! computes the derivative polynomial as the initial h
    ! polynomial and computes l1 no-shift h polynomials.

    implicit none

    integer,intent(in) :: l1

    integer :: i , j , jj , n , nm1
    real(wp) :: xni , t1 , t2

    n = nn - 1
    nm1 = n - 1
    do i = 1 , n
       xni = nn - i
       hr(i) = xni*pr(i)/real(n,wp)
       hi(i) = xni*pi(i)/real(n,wp)
    enddo
    do jj = 1 , l1
       if ( cmod(hr(n),hi(n))<=eta*10.0_wp*cmod(pr(n),pi(n)) ) then
          ! if the constant term is essentially zero, shift h coefficients.
          do i = 1 , nm1
             j = nn - i
             hr(j) = hr(j-1)
             hi(j) = hi(j-1)
          enddo
          hr(1) = 0.0_wp
          hi(1) = 0.0_wp
       else
          call cdivid(-pr(nn),-pi(nn),hr(n),hi(n),tr,ti)
          do i = 1 , nm1
             j = nn - i
             t1 = hr(j-1)
             t2 = hi(j-1)
             hr(j) = tr*t1 - ti*t2 + pr(j)
             hi(j) = tr*t2 + ti*t1 + pi(j)
          enddo
          hr(1) = pr(1)
          hi(1) = pi(1)
       endif
    enddo
end subroutine noshft

subroutine fxshft(l2,zr,zi,conv)

    ! computes l2 fixed-shift h polynomials and tests for
    ! convergence.
    ! initiates a variable-shift iteration and returns with the
    ! approximate zero if successful.

    implicit none

    integer,intent(in) :: l2 !! limit of fixed shift steps
    real(wp) :: zr , zi !! approximate zero if conv is .true.
    logical :: conv !! logical indicating convergence of stage 3 iteration

    integer :: i , j , n
    real(wp) :: otr , oti , svsr , svsi
    logical :: test , pasd , bool

    n = nn - 1
    ! evaluate p at s.
    call polyev(nn,sr,si,pr,pi,qpr,qpi,pvr,pvi)
    test = .true.
    pasd = .false.
    ! calculate first t = -p(s)/h(s).
    call calct(bool)
    ! main loop for one second stage step.
    do j = 1 , l2
       otr = tr
       oti = ti
        ! compute next h polynomial and new t.
       call nexth(bool)
       call calct(bool)
       zr = sr + tr
       zi = si + ti
        ! test for convergence unless stage 3 has failed once or this
        ! is the last h polynomial .
       if ( .not.(bool .or. .not.test .or. j==l2) ) then
          if ( cmod(tr-otr,ti-oti)>=0.5_wp*cmod(zr,zi) ) then
             pasd = .false.
          elseif ( .not.pasd ) then
             pasd = .true.
          else
            ! the weak convergence test has been passed twice, start the
            ! third stage iteration, after saving the current h polynomial
            ! and shift.
             do i = 1 , n
                shr(i) = hr(i)
                shi(i) = hi(i)
             enddo
             svsr = sr
             svsi = si
             call vrshft(10,zr,zi,conv)
             if ( conv ) return
            ! the iteration failed to converge. turn off testing and restore
            ! h,s,pv and t.
             test = .false.
             do i = 1 , n
                hr(i) = shr(i)
                hi(i) = shi(i)
             enddo
             sr = svsr
             si = svsi
             call polyev(nn,sr,si,pr,pi,qpr,qpi,pvr,pvi)
             call calct(bool)
          endif
       endif
    enddo
    ! attempt an iteration with final h polynomial from second stage.
    call vrshft(10,zr,zi,conv)

end subroutine fxshft

subroutine vrshft(l3,zr,zi,conv)

    ! carries out the third stage iteration.

    implicit none

    integer,intent(in) :: l3 !! limit of steps in stage 3.
    real(wp) :: zr , zi !! on entry contains the initial iterate, if the
                        !! iteration converges it contains the final iterate
                        !! on exit.
    logical :: conv !! .true. if iteration converges

    real(wp) :: mp , ms , omp , relstp , r1 , r2 , tp
    logical :: b , bool
    integer :: i , j

    conv = .false.
    b = .false.
    sr = zr
    si = zi

    ! main loop for stage three
    do i = 1 , l3
        ! evaluate p at s and test for convergence.
       call polyev(nn,sr,si,pr,pi,qpr,qpi,pvr,pvi)
       mp = cmod(pvr,pvi)
       ms = cmod(sr,si)
       if ( mp>20.0_wp*errev(nn,qpr,qpi,ms,mp,are,mre) ) then
          if ( i==1 ) then
             omp = mp
          elseif ( b .or. mp<omp .or. relstp>=0.05_wp ) then
            ! exit if polynomial value increases significantly.
             if ( mp*0.1_wp>omp ) return
             omp = mp
          else
            ! iteration has stalled. probably a cluster of zeros. do 5 fixed
            ! shift steps into the cluster to force one zero to dominate.
             tp = relstp
             b = .true.
             if ( relstp<eta ) tp = eta
             r1 = sqrt(tp)
             r2 = sr*(1.0_wp+r1) - si*r1
             si = sr*r1 + si*(1.0_wp+r1)
             sr = r2
             call polyev(nn,sr,si,pr,pi,qpr,qpi,pvr,pvi)
             do j = 1 , 5
                call calct(bool)
                call nexth(bool)
             enddo
             omp = infin
          endif
            ! calculate next iterate.
          call calct(bool)
          call nexth(bool)
          call calct(bool)
          if ( .not.(bool) ) then
             relstp = cmod(tr,ti)/cmod(sr,si)
             sr = sr + tr
             si = si + ti
          endif
       else
          ! polynomial value is smaller in value than a bound on the error
          ! in evaluating p, terminate the iteration.
          conv = .true.
          zr = sr
          zi = si
          return
       endif
    enddo

end subroutine vrshft

subroutine calct(bool)

    ! computes `t = -p(s)/h(s)`.

    implicit none

    logical,intent(out) :: bool !! logical, set true if `h(s)` is essentially zero.

    real(wp) :: hvr , hvi
    integer :: n

    n = nn - 1
    ! evaluate h(s).
    call polyev(n,sr,si,hr,hi,qhr,qhi,hvr,hvi)
    bool = cmod(hvr,hvi)<=are*10.0_wp*cmod(hr(n),hi(n))
    if ( bool ) then
       tr = 0.0_wp
       ti = 0.0_wp
    else
        call cdivid(-pvr,-pvi,hvr,hvi,tr,ti)
    end if

end subroutine calct

subroutine nexth(bool)

    ! calculates the next shifted h polynomial.

    implicit none

    logical,intent(in) :: bool !! logical, if .true. `h(s)` is essentially zero

    real(wp) :: t1 , t2
    integer :: j , n , nm1

    n = nn - 1
    nm1 = n - 1
    if ( bool ) then
       ! if h(s) is zero replace h with qh.
       do j = 2 , n
          hr(j) = qhr(j-1)
          hi(j) = qhi(j-1)
       enddo
       hr(1) = 0.0_wp
       hi(1) = 0.0_wp
    else
        do j = 2 , n
            t1 = qhr(j-1)
            t2 = qhi(j-1)
            hr(j) = tr*t1 - ti*t2 + qpr(j)
            hi(j) = tr*t2 + ti*t1 + qpi(j)
        enddo
        hr(1) = qpr(1)
        hi(1) = qpi(1)
    end if

end subroutine nexth

subroutine polyev(nn,sr,si,pr,pi,qr,qi,pvr,pvi)

    ! evaluates a polynomial  p  at  s  by the horner recurrence
    ! placing the partial sums in q and the computed value in pv.

    implicit none

    integer,intent(in) :: nn
    real(wp) :: pr(nn) , pi(nn) , qr(nn) , qi(nn) , sr , si , pvr , pvi

    real(wp) :: t
    integer :: i

    qr(1) = pr(1)
    qi(1) = pi(1)
    pvr = qr(1)
    pvi = qi(1)
    do i = 2 , nn
       t = pvr*sr - pvi*si + pr(i)
       pvi = pvr*si + pvi*sr + pi(i)
       pvr = t
       qr(i) = pvr
       qi(i) = pvi
    enddo

end subroutine polyev

real(wp) function errev(nn,qr,qi,ms,mp,are,mre)

    ! bounds the error in evaluating the polynomial
    ! by the horner recurrence.

    implicit none

    integer,intent(in) :: nn
    real(wp) :: qr(nn), qi(nn) !! the partial sums
    real(wp) :: ms !! modulus of the point
    real(wp) :: mp !! modulus of polynomial value
    real(wp) :: are, mre !! error bounds on complex addition and multiplication

    integer :: i
    real(wp) :: e

    e = cmod(qr(1),qi(1))*mre/(are+mre)
    do i = 1 , nn
       e = e*ms + cmod(qr(i),qi(i))
    enddo
    errev = e*(are+mre) - mp*mre

end function errev

real(wp) function cauchy(nn,pt,q)

    ! cauchy computes a lower bound on the moduli of
    ! the zeros of a polynomial

    implicit none

    integer,intent(in) :: nn
    real(wp) :: q(nn)
    real(wp) :: pt(nn) !! the modulus of the coefficients.

    integer :: i , n
    real(wp) :: x , xm , f , dx , df

    pt(nn) = -pt(nn)
    ! compute upper estimate of bound.
    n = nn - 1
    x = exp((log(-pt(nn))-log(pt(1)))/real(n,wp))
    if ( pt(n)/=0.0_wp ) then
       ! if newton step at the origin is better, use it.
       xm = -pt(nn)/pt(n)
       if ( xm<x ) x = xm
    endif

    do
        ! chop the interval (0,x) unitl f <= 0.
        xm = x*0.1_wp
        f = pt(1)
        do i = 2 , nn
           f = f*xm + pt(i)
        enddo
        if ( f<=0.0_wp ) then
           dx = x
           do
               ! newton iteration until x converges to two decimal places.
               if ( abs(dx/x)<=0.005_wp ) then
                  cauchy = x
                  exit
               end if
               q(1) = pt(1)
               do i = 2 , nn
                  q(i) = q(i-1)*x + pt(i)
               enddo
               f = q(nn)
               df = q(1)
               do i = 2 , n
                  df = df*x + q(i)
               enddo
               dx = f/df
               x = x - dx
           end do
           exit
        else
           x = xm
        endif
    end do

end function cauchy

real(wp) function scale(nn,pt,eta,infin,smalno,base)

    ! returns a scale factor to multiply the coefficients of the
    ! polynomial. the scaling is done to avoid overflow and to avoid
    ! undetected underflow interfering with the convergence
    ! criterion.  the factor is a power of the base.

    implicit none

    integer :: nn
    real(wp) :: pt(nn)  !! modulus of coefficients of p
    real(wp) :: eta , infin , smalno , base !! constants describing the
                                            !! floating point arithmetic.

    real(wp) :: hi , lo , max , min , x , sc
    integer :: i , l

    ! find largest and smallest moduli of coefficients.
    hi = sqrt(infin)
    lo = smalno/eta
    max = 0.0_wp
    min = infin
    do i = 1 , nn
       x = pt(i)
       if ( x>max ) max = x
       if ( x/=0.0_wp .and. x<min ) min = x
    enddo
    ! scale only if there are very large or very small components.
    scale = 1.0_wp
    if ( min>=lo .and. max<=hi ) return
    x = lo/min
    if ( x>1.0_wp ) then
       sc = x
       if ( infin/sc>max ) sc = 1.0_wp
    else
       sc = 1.0_wp/(sqrt(max)*sqrt(min))
    endif
    l = log(sc)/log(base) + 0.5_wp
    scale = base**l

end function scale

subroutine cdivid(ar,ai,br,bi,cr,ci)

    ! complex division c = a/b, avoiding overflow.

    implicit none

    real(wp) :: ar , ai , br , bi , cr , ci , r , d

    if ( br==0.0_wp .and. bi==0.0_wp ) then
       ! division by zero, c = infinity.
       cr = infin
       ci = infin
    elseif ( abs(br)>=abs(bi) ) then
       r = bi/br
       d = br + r*bi
       cr = (ar+ai*r)/d
       ci = (ai-ar*r)/d
    else
        r = br/bi
        d = bi + r*br
        cr = (ar*r+ai)/d
        ci = (ai*r-ar)/d
    end if

end subroutine cdivid

real(wp) function cmod(r,i)

    implicit none

    real(wp) :: r , i , ar , ai

    ar = abs(r)
    ai = abs(i)
    if ( ar<ai ) then
       cmod = ai*sqrt(1.0_wp+(ar/ai)**2)
    elseif ( ar<=ai ) then
       cmod = ar*sqrt(2.0_wp)
    else
       cmod = ar*sqrt(1.0_wp+(ai/ar)**2)
    end if

end function cmod

end subroutine cpoly
!*****************************************************************************************

!*****************************************************************************************
!>
!  Numerical computation of the roots of a polynomial having
!  complex coefficients, based on aberth's method.
!
!  this routine approximates the roots of the polynomial
!  `p(x)=a(n+1)x^n+a(n)x^(n-1)+...+a(1), a(j)=cr(j)+i ci(j), i**2=-1`,
!  where `a(1)` and `a(n+1)` are nonzero.
!
!  the coefficients are complex numbers. the routine is fast, robust
!  against overflow, and allows to deal with polynomials of any degree.
!  overflow situations are very unlikely and may occurr if there exist
!  simultaneously coefficients of moduli close to big and close to
!  small, i.e., the greatest and the smallest positive real(wp) numbers,
!  respectively. in this limit situation the program outputs a warning
!  message. the computation can be speeded up by performing some side
!  computations in single precision, thus slightly reducing the
!  robustness of the program (see the comments in the routine aberth).
!  besides a set of approximations to the roots, the program delivers a
!  set of a-posteriori error bounds which are guaranteed in the most
!  part of cases. in the situation where underflow does not allow to
!  compute a guaranteed bound, the program outputs a warning message
!  and sets the bound to 0. in the situation where the root cannot be
!  represented as a complex(wp) number the error bound is set to -1.
!
!  the computation is performed by means of aberth's method
!  according to the formula
!```
!           x(i)=x(i)-newt/(1-newt*abcorr), i=1,...,n             (1)
!```
!  where `newt=p(x(i))/p'(x(i))` is the newton correction and `abcorr=
!  =1/(x(i)-x(1))+...+1/(x(i)-x(i-1))+1/(x(i)-x(i+1))+...+1/(x(i)-x(n))`
!  is the aberth correction to the newton method.
!
!  the value of the newton correction is computed by means of the
!  synthetic division algorithm (ruffini-horner's rule) if |x|<=1,
!  otherwise the following more robust (with respect to overflow)
!  formula is applied:
!```
!                    newt=1/(n*y-y**2 r'(y)/r(y))                 (2)
!```
!  where
!```
!                    y=1/x
!                    r(y)=a(1)*y**n+...+a(n)*y+a(n+1)            (2')
!```
!  this computation is performed by the routine [[newton]].
!
!  the starting approximations are complex numbers that are
!  equispaced on circles of suitable radii. the radius of each
!  circle, as well as the number of roots on each circle and the
!  number of circles, is determined by applying rouche's theorem
!  to the functions `a(k+1)*x**k` and `p(x)-a(k+1)*x**k, k=0,...,n`.
!  this computation is performed by the routine [[start]].
!
!### stop condition
!
! if the condition
!```
!                     |p(x(j))|<eps s(|x(j)|)                      (3)
!```
! is satisfied, where `s(x)=s(1)+x*s(2)+...+x**n * s(n+1)`,
! `s(i)=|a(i)|*(1+3.8*(i-1))`, `eps` is the machine precision (eps=2**-53
! for the ieee arithmetic), then the approximation `x(j)` is not updated
! and the subsequent iterations (1)  for `i=j` are skipped.
! the program stops if the condition (3) is satisfied for `j=1,...,n`,
! or if the maximum number `nitmax` of iterations has been reached.
! the condition (3) is motivated by a backward rounding error analysis
! of the ruffini-horner rule, moreover the condition (3) guarantees
! that the computed approximation `x(j)` is an exact root of a slightly
! perturbed polynomial.
!
!### inclusion disks, a-posteriori error bounds
!
! for each approximation `x` of a root, an a-posteriori absolute error
! bound r is computed according to the formula
!```
!                   r=n(|p(x)|+eps s(|x|))/|p'(x)|                 (4)
!```
! this provides an inclusion disk of center `x` and radius `r` containing a
! root.
!
!### Reference
!  * Dario Andrea Bini, "[Numerical computation of polynomial zeros by means of Aberth's method](https://link.springer.com/article/10.1007/BF02207694)"
!    Numerical Algorithms volume 13, pages 179-200 (1996)
!
!### History
!  * version 1.4, june 1996
!    (d. bini, dipartimento di matematica, universita' di pisa)
!    (bini@dm.unipi.it)
!    work performed under the support of the esprit bra project 6846 posso
!    Source: [Netlib](https://netlib.org/numeralgo/na10)
!  * Jacob Williams, 9/19/2022, modernized this code

    subroutine polzeros(n, poly, nitmax, root, radius, err)

        implicit none

        integer,intent(in) :: n !! degree of the polynomial.
        complex(wp),intent(in) :: poly(n + 1) !! complex vector of n+1 components, `poly(i)` is the
                                              !! coefficient of `x**(i-1), i=1,...,n+1` of the polynomial `p(x)`
        integer,intent(in) :: nitmax !! the max number of allowed iterations.
        complex(wp),intent(out) :: root(n) !! complex vector of `n` components, containing the
                                           !! approximations to the roots of `p(x)`.
        real(wp),intent(out) :: radius(n) !! real vector of `n` components, containing the error bounds to
                                          !! the approximations of the roots, i.e. the disk of center
                                          !! `root(i)` and radius `radius(i)` contains a root of `p(x)`, for
                                          !! `i=1,...,n`. `radius(i)` is set to -1 if the corresponding root
                                          !! cannot be represented as floating point due to overflow or
                                          !! underflow.
        logical,intent(out) :: err(n)     !! vector of `n` components detecting an error condition:
                                          !!
                                          !!  * `err(j)=.true.` if after `nitmax` iterations the stop condition
                                          !!    (3) is not satisfied for x(j)=root(j);
                                          !!  * `err(j)=.false.`  otherwise, i.e., the root is reliable,
                                          !!    i.e., it can be viewed as an exact root of a
                                          !!    slightly perturbed polynomial.
                                          !!
                                          !! the vector `err` is used also in the routine convex hull for
                                          !! storing the abscissae of the vertices of the convex hull.

        integer :: iter !! number of iterations peformed
        real(wp) :: apoly(n + 1) !! auxiliary variable: real vector of n+1 components used to store the moduli of
                                 !! the coefficients of p(x) and the coefficients of s(x) used
                                 !! to test the stop condition (3).
        real(wp) :: apolyr(n + 1) !! auxiliary variable: real vector of n+1 components used to test the stop
                                  !! condition
        integer :: i, nzeros
        complex(wp) :: corr, abcorr
        real(wp) :: amax

        real(wp),parameter :: eps   = epsilon(1.0_wp)
        real(wp),parameter :: small = tiny(1.0_wp)
        real(wp),parameter :: big   = huge(1.0_wp)

        ! check consistency of data
        if (abs(poly(n + 1)) == 0.0_wp) then
            error stop 'inconsistent data: the leading coefficient is zero'
        end if
        if (abs(poly(1)) == 0.0_wp) then
            error stop 'the constant term is zero: deflate the polynomial'
        end if
        ! compute the moduli of the coefficients
        amax = 0.0_wp
        do i = 1, n + 1
            apoly(i) = abs(poly(i))
            amax = max(amax, apoly(i))
            apolyr(i) = apoly(i)
        end do
        if ((amax) >= (big/(n + 1))) then
            write (*, *) 'warning: coefficients too big, overflow is likely'
        end if
        ! initialize
        do i = 1, n
            radius(i) = 0.0_wp
            err(i) = .true.
        end do
        ! select the starting points
        call start(n, apolyr, root, radius, nzeros, small, big)
        ! compute the coefficients of the backward-error polynomial
        do i = 1, n + 1
            apolyr(n - i + 2) = eps*apoly(i)*(3.8_wp*(n - i + 1) + 1)
            apoly(i) = eps*apoly(i)*(3.8_wp*(i - 1) + 1)
        end do
        if ((apoly(1) == 0.0_wp) .or. (apoly(n + 1) == 0.0_wp)) then
            write (*, *) 'warning: the computation of some inclusion radius'
            write (*, *) 'may fail. this is reported by radius=0'
        end if
        do i = 1, n
            err(i) = .true.
            if (radius(i) == -1) err(i) = .false.
        end do
        ! starts aberth's iterations
        do iter = 1, nitmax
          do i = 1, n
              if (err(i)) then
                  call newton(n, poly, apoly, apolyr, root(i), small, radius(i), corr, err(i))
                  if (err(i)) then
                      call aberth(n, i, root, abcorr)
                      root(i) = root(i) - corr/(1 - corr*abcorr)
                  else
                      nzeros = nzeros + 1
                      if (nzeros == n) return
                  end if
              end if
          end do
        end do

    end subroutine polzeros

    subroutine newton(n, poly, apoly, apolyr, z, small, radius, corr, again)

        !! compute the newton's correction, the inclusion radius (4) and checks
        !! the stop condition (3)

        implicit none

        integer,intent(in) :: n !! degree of the polynomial p(x)
        complex(wp),intent(in) :: poly(n + 1) !! coefficients of the polynomial p(x)
        real(wp),intent(in) :: apoly(n + 1) !! upper bounds on the backward perturbations on the
                                            !! coefficients of p(x) when applying ruffini-horner's rule
        real(wp),intent(in) :: apolyr(n + 1) !! upper bounds on the backward perturbations on the
                                             !! coefficients of p(x) when applying (2), (2')
        complex(wp),intent(in) :: z !! value at which the newton correction is computed
        real(wp),intent(in) :: small !! the min positive real(wp), small=2**(-1074) for the ieee.
        real(wp),intent(out) :: radius !! upper bound to the distance of z from the closest root of
                                       !! the polynomial computed according to (4).
        complex(wp),intent(out) :: corr !! newton's correction
        logical,intent(out) :: again !! this variable is .true. if the computed value p(z) is
                                     !! reliable, i.e., (3) is not satisfied in z. again is
                                     !! .false., otherwise.

        integer :: i
        complex(wp) :: p, p1, zi, den, ppsp
        real(wp) :: ap, az, azi, absp

        az = abs(z)
        ! if |z|<=1 then apply ruffini-horner's rule for p(z)/p'(z)
        ! and for the computation of the inclusion radius
        if (az <= 1.0_wp) then
            p = poly(n + 1)
            ap = apoly(n + 1)
            p1 = p
            do i = n, 2, -1
                p = p*z + poly(i)
                p1 = p1*z + p
                ap = ap*az + apoly(i)
            end do
            p = p*z + poly(1)
            ap = ap*az + apoly(1)
            corr = p/p1
            absp = abs(p)
            ap = ap
            again = (absp > (small + ap))
            if (.not. again) radius = n*(absp + ap)/abs(p1)
        else
            ! if |z|>1 then apply ruffini-horner's rule to the reversed polynomial
            ! and use formula (2) for p(z)/p'(z). analogously do for the inclusion
            ! radius.
            zi = 1.0_wp/z
            azi = 1.0_wp/az
            p = poly(1)
            p1 = p
            ap = apolyr(n + 1)
            do i = n, 2, -1
                p = p*zi + poly(n - i + 2)
                p1 = p1*zi + p
                ap = ap*azi + apolyr(i)
            end do
            p = p*zi + poly(n + 1)
            ap = ap*azi + apolyr(1)
            absp = abs(p)
            again = (absp > (small + ap))
            ppsp = (p*z)/p1
            den = n*ppsp - 1
            corr = z*(ppsp/den)
            if (again) return
            radius = abs(ppsp) + (ap*az)/abs(p1)
            radius = n*radius/abs(den)
            radius = radius*az
        end if

    end subroutine newton

    subroutine aberth(n, j, root, abcorr)

        !! compute the aberth correction. to save time, the reciprocation of
        !! root(j)-root(i) could be performed in single precision (complex*8)
        !! in principle this might cause overflow if both root(j) and root(i)
        !! have too small moduli.

        implicit none

        integer,intent(in) :: n !! degree of the polynomial
        integer,intent(in) :: j  !! index of the component of root with respect to which the
                                 !! aberth correction is computed
        complex(wp),intent(in) :: root(n) !! vector containing the current approximations to the roots
        complex(wp),intent(out) :: abcorr !! aberth's correction (compare (1))

        integer :: i
        complex(wp) :: z, zj

        abcorr = 0.0_wp
        zj = root(j)
        do i = 1, j - 1
            z = zj - root(i)
            abcorr = abcorr + 1.0_wp/z
        end do
        do i = j + 1, n
            z = zj - root(i)
            abcorr = abcorr + 1.0_wp/z
        end do

    end subroutine aberth

    subroutine start(n, a, y, radius, nz, small, big)

        !! compute the starting approximations of the roots
        !!
        !! this routines selects starting approximations along circles center at
        !! 0 and having suitable radii. the computation of the number of circles
        !! and of the corresponding radii is performed by computing the upper
        !! convex hull of the set (i,log(a(i))), i=1,...,n+1.

        implicit none

        integer,intent(in) :: n !! number of the coefficients of the polynomial
        real(wp),intent(inout) :: a(n + 1) !! moduli of the coefficients of the polynomial
        complex(wp),intent(out) :: y(n) !! starting approximations
        real(wp),intent(out) :: radius(n) !! if a component is -1 then the corresponding root has a
                                          !! too big or too small modulus in order to be represented
                                          !! as double float with no overflow/underflow
        integer,intent(out) :: nz !! number of roots which cannot be represented without
                                  !! overflow/underflow
        real(wp),intent(in) :: small !! the min positive real(wp), small=2**(-1074) for the ieee.
        real(wp),intent(in) :: big !! the max real(wp), big=2**1023 for the ieee standard.

        logical :: h(n + 1) !! auxiliary variable: needed for the computation of the convex hull

        integer :: i, iold, nzeros, j, jj
        real(wp) :: r, th, ang, temp
        real(wp) :: xsmall, xbig

        real(wp),parameter :: pi2 = 2.0_wp * pi
        real(wp),parameter :: sigma = 0.7_wp

        xsmall = log(small)
        xbig = log(big)
        nz = 0
        ! compute the logarithm a(i) of the moduli of the coefficients of
        ! the polynomial and then the upper covex hull of the set (a(i),i)
        do i = 1, n + 1
            if (a(i) /= 0.0_wp) then
                a(i) = log(a(i))
            else
                a(i) = -1.0e30_wp ! maybe replace with -huge(1.0_wp) ?? -JW
            end if
        end do
        call cnvex(n + 1, a, h)
        ! given the upper convex hull of the set (a(i),i) compute the moduli
        ! of the starting approximations by means of rouche's theorem
        iold = 1
        th = pi2/n
        do i = 2, n + 1
            if (h(i)) then
                nzeros = i - iold
                temp = (a(iold) - a(i))/nzeros
                ! check if the modulus is too small
                if ((temp < -xbig) .and. (temp >= xsmall)) then
                    write (*, *) 'warning:', nzeros, ' zero(s) are too small to'
                    write (*, *) 'represent their inverses as complex(wp), they'
                    write (*, *) 'are replaced by small numbers, the corresponding'
                    write (*, *) 'radii are set to -1'
                    nz = nz + nzeros
                    r = 1.0_wp/big
                end if
                if (temp < xsmall) then
                    nz = nz + nzeros
                    write (*, *) 'warning: ', nzeros, ' zero(s) are too small to be'
                    write (*, *) 'represented as complex(wp), they are set to 0'
                    write (*, *) 'the corresponding radii are set to -1'
                end if
                ! check if the modulus is too big
                if (temp > xbig) then
                    r = big
                    nz = nz + nzeros
                    write (*, *) 'warning: ', nzeros, ' zeros(s) are too big to be'
                    write (*, *) 'represented as complex(wp),'
                    write (*, *) 'the corresponding radii are set to -1'
                end if
                if ((temp <= xbig) .and. (temp > max(-xbig, xsmall))) r = exp(temp)
                ! compute nzeros approximations equally distributed in the disk of
                ! radius r
                ang = pi2/nzeros
                do j = iold, i - 1
                    jj = j - iold + 1
                    if ((r <= (1.0_wp/big)) .or. (r == big)) radius(j) = -1.0_wp
                    y(j) = r*(cos(ang*jj + th*i + sigma) + cmplx(0.0_wp, 1.0_wp, wp)*sin(ang*jj + th*i + sigma))
                end do
                iold = i
            end if
        end do

    end subroutine start

    subroutine cnvex(n, a, h)

        !! compute the upper convex hull of the set (i,a(i)), i.e., the set of
        !! vertices (i_k,a(i_k)), k=1,2,...,m, such that the points (i,a(i)) lie
        !! below the straight lines passing through two consecutive vertices.
        !! the abscissae of the vertices of the convex hull equal the indices of
        !! the true  components of the logical output vector h.
        !! the used method requires o(nlog n) comparisons and is based on a
        !! divide-and-conquer technique. once the upper convex hull of two
        !! contiguous sets  (say, {(1,a(1)),(2,a(2)),...,(k,a(k))} and
        !! {(k,a(k)), (k+1,a(k+1)),...,(q,a(q))}) have been computed, then
        !! the upper convex hull of their union is provided by the subroutine
        !! cmerge. the program starts with sets made up by two consecutive
        !! points, which trivially constitute a convex hull, then obtains sets
        !! of 3,5,9... points,  up to  arrive at the entire set.
        !! the program uses the subroutine  cmerge; the subroutine cmerge uses
        !! the subroutines left, right and ctest. the latter tests the convexity
        !! of the angle formed by the points (i,a(i)), (j,a(j)), (k,a(k)) in the
        !! vertex (j,a(j)) up to within a given tolerance toler, where i<j<k.

        implicit none

        integer,intent(in) :: n
        real(wp) :: a(n)
        logical,intent(out) :: h(n)

        integer :: i, j, k, m, nj, jc

        do i = 1, n
            h(i) = .true.
        end do
        ! compute k such that n-2<=2**k<n-1
        k = int(log(n - 2.0_wp)/log(2.0_wp))
        if (2**(k + 1) <= (n - 2)) k = k + 1
        ! for each m=1,2,4,8,...,2**k, consider the nj pairs of consecutive
        ! sets made up by m+1 points having the common vertex
        ! (jc,a(jc)), where jc=m*(2*j+1)+1 and j=0,...,nj,
        ! nj=max(0,int((n-2-m)/(m+m))).
        ! compute the upper convex hull of their union by means of the
        ! subroutine cmerge
        m = 1
        do i = 0, k
            nj = max(0, int((n - 2 - m)/(m + m)))
            do j = 0, nj
                jc = (j + j + 1)*m + 1
                call cmerge(n, a, jc, m, h)
            end do
            m = m + m
        end do

    end subroutine cnvex

    subroutine left(n, h, i, il)

        !! given as input the integer i and the vector h of logical, compute the
        !! the maximum integer il such that il<i and h(il) is true.

        implicit none

        integer,intent(in) :: n !! length of the vector h
        integer,intent(in) :: i !! integer
        logical,intent(in) :: h(n) !! vector of logical
        integer,intent(out) :: il !! maximum integer such that il<i, h(il)=.true.

        do il = i - 1, 0, -1
            if (h(il)) return
        end do

    end subroutine left

    subroutine right(n, h, i, ir)

        !! given as input the integer i and the vector h of logical, compute the
        !! the minimum integer ir such that ir>i and h(il) is true.

        implicit none

        integer,intent(in) :: n !! length of the vector h
        logical ,intent(in):: h(n) !! vector of logical
        integer,intent(in) :: i !! integer
        integer,intent(out) :: ir !! minimum integer such that ir>i, h(ir)=.true.

        do ir = i + 1, n
            if (h(ir)) return
        end do

    end subroutine right

    subroutine cmerge(n, a, i, m, h)

        !! given the upper convex hulls of two consecutive sets of pairs
        !! (j,a(j)), compute the upper convex hull of their union

        implicit none

        integer,intent(in) :: n !! length of the vector a
        real(wp),intent(in) :: a(n) !! vector defining the points (j,a(j))
        integer,intent(in) :: i !! abscissa of the common vertex of the two sets
        integer,intent(in) :: m !! the number of elements of each set is m+1
        logical,intent(out) :: h(n) !! vector defining the vertices of the convex hull, i.e.,
                                    !! h(j) is .true. if (j,a(j)) is a vertex of the convex hull
                                    !! this vector is used also as output.

        integer :: ir, il, irr, ill
        logical :: tstl, tstr

        ! at the left and the right of the common vertex (i,a(i)) determine
        ! the abscissae il,ir, of the closest vertices of the upper convex
        ! hull of the left and right sets, respectively
        call left(n, h, i, il)
        call right(n, h, i, ir)
        ! check the convexity of the angle formed by il,i,ir
        if (ctest(n, a, il, i, ir)) then
            return
        else
            ! continue the search of a pair of vertices in the left and right
            ! sets which yield the upper convex hull
            h(i) = .false.
            do
                if (il == (i - m)) then
                    tstl = .true.
                else
                    call left(n, h, il, ill)
                    tstl = ctest(n, a, ill, il, ir)
                end if
                if (ir == min(n, i + m)) then
                    tstr = .true.
                else
                    call right(n, h, ir, irr)
                    tstr = ctest(n, a, il, ir, irr)
                end if
                h(il) = tstl
                h(ir) = tstr
                if (tstl .and. tstr) return
                if (.not. tstl) il = ill
                if (.not. tstr) ir = irr
            end do
        end if

    end subroutine cmerge

    function ctest(n, a, il, i, ir)

        !! test the convexity of the angle formed by (il,a(il)), (i,a(i)),
        !! (ir,a(ir)) at the vertex (i,a(i)), up to within the tolerance
        !! toler. if convexity holds then the function is set to .true.,
        !! otherwise ctest=.false. the parameter toler is set to 0.4 by default.

        implicit none

        integer,intent(in) :: n !! length of the vector a
        integer,intent(in) :: i !! integers such that il<i<ir
        integer,intent(in) :: il !! integers such that il<i<ir
        integer,intent(in) :: ir !! integers such that il<i<ir
        real(wp),intent(in) :: a(n) !! vector of double
        logical :: ctest !! Result:
                         !!
                         !! * .true. if the angle formed by (il,a(il)), (i,a(i)), (ir,a(ir)) at
                         !!   the vertex (i,a(i)), is convex up to within the tolerance
                         !!   toler, i.e., if
                         !!   (a(i)-a(il))*(ir-i)-(a(ir)-a(i))*(i-il)>toler.
                         !!
                         !! * .false.,  otherwise.

        real(wp) :: s1, s2

        real(wp), parameter :: toler = 0.4_wp

        s1 = a(i) - a(il)
        s2 = a(ir) - a(i)
        s1 = s1*(ir - i)
        s2 = s2*(i - il)
        ctest = .false.
        if (s1 > (s2 + toler)) ctest = .true.

    end function ctest
!*****************************************************************************************

!*****************************************************************************************
!>
!  FPML: Fourth order Parallelizable Modification of Laguerre's method
!
!### Reference
!  * Thomas R. Cameron,
!    "An effective implementation of a modified Laguerre method for the roots of a polynomial",
!    Numerical Algorithms volume 82, pages 1065-1084 (2019)
!    [link](https://link.springer.com/article/10.1007/s11075-018-0641-9)
!
!### History
!  * Author: Thomas R. Cameron, Davidson College
!    Last Modified: 1 November 2018
!    Original code: https://github.com/trcameron/FPML
!  * Jacob Williams, 9/21/2022 : refactored this code a bit.

    subroutine fpml(poly, deg, roots, berr, cond, conv, itmax)

        implicit none

        integer, intent(in)      :: deg !! polynomial degree
        complex(wp), intent(in)  :: poly(deg+1) !! coefficients
        complex(wp), intent(out) :: roots(:) !! the root approximations
        real(wp), intent(out)    :: berr(:) !! the backward error in each approximation
        real(wp), intent(out)    :: cond(:) !! the condition number of each root approximation
        integer, intent(out)     :: conv(:)
        integer, intent(in)      :: itmax

        integer                    :: i, j, nz
        real(wp)                   :: r
        real(wp), dimension(deg+1) :: alpha
        complex(wp)                :: b, c, z

        real(wp), parameter :: big = huge(1.0_wp)
        real(wp), parameter :: small = tiny(1.0_wp)

        ! precheck
        alpha = abs(poly)
        if (alpha(deg+1)<small) then
            write(*,'(A)') 'Warning: leading coefficient too small.'
            return
        elseif (deg==1) then
            roots(1) = -poly(1)/poly(2)
            conv = 1
            berr = 0.0_wp
            cond(1) = (alpha(1) + alpha(2)*abs(roots(1)))/(abs(roots(1))*alpha(2))
            return
        elseif (deg==2) then
            b = -poly(2)/(2.0_wp*poly(3))
            c = sqrt(poly(2)**2-4.0_wp*poly(3)*poly(1))/(2.0_wp*poly(3))
            roots(1) = b - c
            roots(2) = b + c
            conv = 1
            berr = 0.0_wp
            cond(1) = (alpha(1)+alpha(2)*abs(roots(1))+&
                       alpha(3)*abs(roots(1))**2)/(abs(roots(1))*&
                       abs(poly(2)+2.0_wp*poly(3)*roots(1)))
            cond(2) = (alpha(1)+alpha(2)*abs(roots(2))+&
                       alpha(3)*abs(roots(2))**2)/(abs(roots(2))*&
                       abs(poly(2)+2.0_wp*poly(3)*roots(2)))
            return
        end if
        ! initial estimates
        conv = [(0, i=1,deg)]
        nz = 0
        call estimates(alpha, deg, roots, conv, nz)
        ! main loop
        alpha = [(alpha(i)*(3.8_wp*(i-1)+1),i=1,deg+1)]
        main : do i=1,itmax
            do j=1,deg
                if (conv(j)==0) then
                    z = roots(j)
                    r = abs(z)
                    if (r > 1.0_wp) then
                        call rcheck_lag(poly, alpha, deg, b, c, z, r, conv(j), berr(j), cond(j))
                    else
                        call check_lag(poly, alpha, deg, b, c, z, r, conv(j), berr(j), cond(j))
                    end if
                    if (conv(j)==0) then
                        call modify_lag(deg, b, c, z, j, roots)
                        roots(j) = roots(j) - c
                    else
                        nz = nz + 1
                        if (nz==deg) exit main
                    end if
                end if
            end do
        end do main
        ! final check
        if (minval(conv)==1) then
            return
        else
            ! display warrning
            write(*,'(A)') 'Some root approximations did not converge or experienced overflow/underflow.'
            ! compute backward error and condition number for roots that did not converge;
            ! note that this may produce overflow/underflow.
            do j=1,deg
                if (conv(j) /= 1) then
                    z = roots(j)
                    r = abs(z)
                    if (r>1.0_wp) then
                        z = 1.0_wp/z
                        r = 1.0_wp/r
                        c = 0.0_wp
                        b = poly(1)
                        berr(j) = alpha(1)
                        do i=2,deg+1
                            c = z*c + b
                            b = z*b + poly(i)
                            berr(j) = r*berr(j) + alpha(i)
                        end do
                        cond(j) = berr(j)/abs(deg*b-z*c)
                        berr(j) = abs(b)/berr(j)
                    else
                        c = 0
                        b = poly(deg+1)
                        berr(j) = alpha(deg+1)
                        do i=deg,1,-1
                            c = z*c + b
                            b = z*b + poly(i)
                            berr(j) = r*berr(j) + alpha(i)
                        end do
                        cond(j) = berr(j)/(r*abs(c))
                        berr(j) = abs(b)/berr(j)
                    end if
                end if
            end do
        end if

    contains

    !************************************************
    ! Computes backward error of root approximation
    ! with moduli greater than 1.
    ! If the backward error is less than eps, then
    ! both backward error and condition number are
    ! computed. Otherwise, the Laguerre correction terms
    ! are computed and stored in variables b and c.
    !************************************************
    subroutine rcheck_lag(p, alpha, deg, b, c, z, r, conv, berr, cond)
        implicit none
        ! argument variables
        integer, intent(in)        :: deg
        integer, intent(out)       :: conv
        real(wp), intent(in)       :: alpha(:), r
        real(wp), intent(out)      :: berr, cond
        complex(wp), intent(in)    :: p(:), z
        complex(wp), intent(out)   :: b, c
        ! local variables
        integer     :: k
        real(wp)    :: rr
        complex(wp) :: a, zz

        ! evaluate polynomial and derivatives
        zz = 1.0_wp/z
        rr = 1.0_wp/r
        a = p(1)
        b = 0
        c = 0
        berr = alpha(1)
        do k=2,deg+1
            c = zz*c + b
            b = zz*b + a
            a = zz*a + p(k)
            berr = rr*berr + alpha(k)
        end do
        ! laguerre correction/ backward error and condition
        if (abs(a)>eps*berr) then
            b = b/a
            c = 2.0_wp*(c/a)
            c = zz**2*(deg-2*zz*b+zz**2*(b**2-c))
            b = zz*(deg-zz*b)
            if (check_nan_inf(b) .or. check_nan_inf(c)) conv = -1
        else
            cond = berr/abs(deg*a-zz*b)
            berr = abs(a)/berr
            conv = 1
        end if
    end subroutine rcheck_lag

    !************************************************
    ! Computes backward error of root approximation
    ! with moduli less than or equal to 1.
    ! If the backward error is less than eps, then
    ! both backward error and condition number are
    ! computed. Otherwise, the Laguerre correction terms
    ! Gj and Hj are computed and stored in variables
    ! b and c, respectively.
    !************************************************
    subroutine check_lag(p, alpha, deg, b, c, z, r, conv, berr, cond)
        implicit none
        ! argument variables
        integer, intent(in)        :: deg
        integer, intent(out)       :: conv
        real(wp), intent(in)       :: alpha(:), r
        real(wp), intent(out)      :: berr, cond
        complex(wp), intent(in)    :: p(:), z
        complex(wp), intent(out)   :: b, c
        ! local variables
        integer                    :: k
        complex(wp)                :: a

        ! evaluate polynomial and derivatives
        a = p(deg+1)
        b = 0.0_wp
        c = 0.0_wp
        berr = alpha(deg+1)
        do k=deg,1,-1
            c = z*c + b
            b = z*b + a
            a = z*a + p(k)
            berr = r*berr + alpha(k)
        end do
        ! laguerre correction/ backward error and condition
        if (abs(a)>eps*berr) then
            b = b/a
            c = b**2 - 2.0_wp*(c/a)
            if (check_nan_inf(b) .or. check_nan_inf(c)) conv = -1
        else
            cond = berr/(r*abs(b))
            berr = abs(a)/berr
            conv = 1
        end if
    end subroutine check_lag

    !************************************************
    ! Computes modified Laguerre correction term of
    ! the jth rooot approximation.
    ! The coefficients of the polynomial of degree
    ! deg are stored in p, all root approximations
    ! are stored in roots. The values b, and c come
    ! from rcheck_lag or check_lag, c will be used
    ! to return the correction term.
    !************************************************
    subroutine modify_lag(deg, b, c, z, j, roots)
        implicit none
        ! argument variables
        integer, intent(in)        :: deg, j
        complex(wp), intent(in)    :: roots(:), z
        complex(wp), intent(inout) :: b, c
        ! local variables
        integer                    :: k
        complex(wp)                :: t

        ! Aberth correction terms
        do k=1,j-1
            t = 1.0_wp/(z - roots(k))
            b = b - t
            c = c - t**2
        end do
        do k=j+1,deg
            t = 1.0_wp/(z - roots(k))
            b = b - t
            c = c - t**2
        end do
        ! Laguerre correction
        t = sqrt((deg-1)*(deg*c-b**2))
        c = b + t
        b = b - t
        if (abs(b)>abs(c)) then
            c = deg/b
        else
            c = deg/c
        end if
    end subroutine modify_lag

    !************************************************
    ! Computes initial estimates for the roots of an
    ! univariate polynomial of degree deg, whose
    ! coefficients moduli are stored in alpha. The
    ! estimates are returned in the array roots.
    ! The computation is performed as follows: First
    ! the set (i,log(alpha(i))) is formed and the
    ! upper envelope of the convex hull of this set
    ! is computed, its indices are returned in the
    ! array h (in descending order). For i=c-1,1,-1
    ! there are h(i) - h(i+1) zeros placed on a
    ! circle of radius alpha(h(i+1))/alpha(h(i))
    ! raised to the 1/(h(i)-h(i+1)) power.
    !************************************************
    subroutine estimates(alpha, deg, roots, conv, nz)
        implicit none

        real(wp), intent(in)       :: alpha(:)
        integer, intent(in)        :: deg
        complex(wp), intent(inout) :: roots(:)
        integer, intent(inout)     :: conv(:)
        integer, intent(inout)     :: nz

        integer                    :: c, i, j, k, nzeros
        real(wp)                   :: a1, a2, ang, r, th
        integer, dimension(deg+1)  :: h
        real(wp), dimension(deg+1) :: a

        real(wp), parameter :: pi2 = 2.0_wp * pi
        real(wp), parameter :: sigma = 0.7_wp

        ! Log of absolute value of coefficients
        do i=1,deg+1
            if (alpha(i)>0) then
                a(i) = log(alpha(i))
            else
                a(i) = -1.0e30_wp
            end if
        end do
        call conv_hull(deg+1, a, h, c)
        k=0
        th=pi2/deg
        ! Initial Estimates
        do i=c-1,1,-1
            nzeros = h(i)-h(i+1)
            a1 = alpha(h(i+1))**(1.0_wp/nzeros)
            a2 = alpha(h(i))**(1.0_wp/nzeros)
            if (a1 <= a2*small) then
                ! r is too small
                r = 0.0_wp
                nz = nz + nzeros
                conv(k+1:k+nzeros) = -1
                roots(k+1:k+nzeros) = cmplx(0.0_wp,0.0_wp,wp)
            else if (a1 >= a2*big) then
                ! r is too big
                r = big
                nz = nz+nzeros
                conv(k+1:k+nzeros) = -1
                ang = pi2/nzeros
                do j=1,nzeros
                    roots(k+j) = r*cmplx(cos(ang*j+th*h(i)+sigma),sin(ang*j+th*h(i)+sigma),wp)
                end do
            else
                ! r is just right
                r = a1/a2
                ang = pi2/nzeros
                do j=1,nzeros
                    roots(k+j) = r*cmplx(cos(ang*j+th*h(i)+sigma),sin(ang*j+th*h(i)+sigma),wp)
                end do
            end if
            k = k+nzeros
        end do
    end subroutine estimates

    !************************************************
    ! Computex upper envelope of the convex hull of
    ! the points in the array a, which has size n.
    ! The number of vertices in the hull is equal to
    ! c, and they are returned in the first c entries
    ! of the array h.
    ! The computation follows Andrew's monotone chain
    ! algorithm: Each consecutive three pairs are
    ! tested via cross to determine if they form
    ! a clockwise angle, if so that current point
    ! is rejected from the returned set.
    !
    !@note The original version of this code had a bug
    !************************************************
    subroutine conv_hull(n, a, h, c)
        implicit none
        ! argument variables
        integer, intent(in)    :: n
        integer, intent(inout) :: c
        integer, intent(inout) :: h(:)
        real(wp), intent(in)   :: a(:)
        ! local variables
        integer :: i

        ! covex hull
        c=0
        do i=n,1,-1
            !do while(c>=2 .and. cross(h, a, c, i)<eps) ! bug in original code here
            do while (c>=2) ! bug in original code here
               if (cross(h, a, c, i)>=eps) exit
                c = c - 1
            end do
            c = c + 1
            h(c) = i
        end do
    end subroutine conv_hull

    !************************************************
    ! Returns 2D cross product of OA and OB vectors,
    ! where
    ! O=(h(c-1),a(h(c-1))),
    ! A=(h(c),a(h(c))),
    ! B=(i,a(i)).
    ! If det>0, then OAB makes counter-clockwise turn.
    !************************************************
    function cross(h, a, c, i) result(det)
        implicit none
        ! argument variables
        integer, intent(in)  :: c, i
        integer, intent(in)  :: h(:)
        real(wp), intent(in) :: a(:)
        ! local variables
        real(wp) :: det

        ! determinant
        det = (a(i)-a(h(c-1)))*(h(c)-h(c-1)) - (a(h(c))-a(h(c-1)))*(i-h(c-1))

    end function cross

    !************************************************
    ! Check if real or imaginary part of complex
    ! number a is either NaN or Inf.
    !************************************************
    function check_nan_inf(a) result(res)
        implicit none
        ! argument variables
        complex(wp),intent(in) :: a
        ! local variables
        logical  :: res
        real(wp) :: re_a, im_a

        ! check for nan and inf
        re_a = real(a,wp)
        im_a = aimag(a)
        res = ieee_is_nan(re_a) .or. ieee_is_nan(im_a) .or. &
              (abs(re_a)>big) .or. (abs(im_a)>big)

    end function check_nan_inf

    end subroutine fpml
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the roots of a cubic equation with real coefficients.
!
!### Reference
!  * V. I. Lebedev, "On formulae for roots of cubic equation",
!    Sov. J. Numer. Anal. Math. Modelling, Vol.6, No.4, pp. 315-324 (1991)
!
!### History
!  * Jacob Williams, 9/23/2022 : based on the `TC` routine in the reference.

    subroutine rroots_chebyshev_cubic(coeffs,zr,zi)

        implicit none

        real(wp),dimension(4),intent(in) :: coeffs !! vector of coefficients in order of decreasing powers
        real(wp), dimension(3), intent(out) :: zr !! output vector of real parts of the zeros
        real(wp), dimension(3), intent(out) :: zi !! output vector of imaginary parts of the zeros

        integer :: l !! number of complex roots (0 or 2)
        real(wp) :: a,b,c,d,p,t1,t2,t3,t4,t,x1,x2,x3

        real(wp),parameter :: sqrt3 = sqrt(3.0_wp)
        real(wp),parameter :: s = 1.0_wp / 3.0_wp
        real(wp),parameter :: small = 10.0_wp**int(log(epsilon(1.0_wp)))  ! this was 1.0d-32 in the original code

        ! coefficients:
        a = coeffs(1)
        b = coeffs(2)
        c = coeffs(3)
        d = coeffs(4)

        main : block
            t = sqrt3
            t2 = b*b
            t3 = 3.0_wp*a
            t4 = t3*c
            p = t2 - t4
            x3 = abs(p)
            x3 = sqrt(x3)
            x1 = b*(t4-p-p) - 3.0_wp*t3*t3*d
            x2 = abs(x1)
            x2 = x2**s
            t2 = 1.0_wp/t3
            t3 = b*t2
            if ( x3>small*x2 ) then
                t1 = 0.5_wp*x1/(p*x3)
                x2 = abs(t1)
                t2 = x3*t2
                t = t*t2
                t4 = x2*x2
                if ( p<0.0_wp ) then
                    p = x2 + sqrt(t4+1.0_wp)
                    p = p**s
                    t4 = -1.0_wp/p
                    if ( t1>=0.0_wp ) t2 = -t2
                    x1 = (p+t4)*t2
                    x2 = -0.5_wp*x1
                    x3 = 0.5_wp*t*(p-t4)
                    l = 2
                    exit main
                else
                    x3 = abs(1.0_wp-t4)
                    x3 = sqrt(x3)
                    if ( t4>1.0_wp ) then
                        p = (x2+x3)**s
                        t4 = 1.0_wp/p
                        if ( t1<0 ) t2 = -t2
                        x1 = (p+t4)*t2
                        x2 = -0.5_wp*x1
                        x3 = 0.5_wp*t*(p-t4)
                        l = 2
                        exit main
                    else
                        t4 = atan2(x3,t1)*s
                        x3 = cos(t4)
                        t4 = sqrt(1.0_wp-x3*x3)*t
                        x3 = x3*t2
                        x1 = x3 + x3
                        x2 = t4 - x3
                        x3 = -(t4+x3)
                        if ( x2<=x3 ) then
                            t2 = x2
                            x2 = x3
                            x3 = t2
                        endif
                    endif
                endif
            else
                if ( x1<0.0_wp ) x2 = -x2
                x1 = x2*t2
                x2 = -0.5_wp*x1
                x3 = -t*x2
                if ( abs(x3)>small ) then
                    l = 2
                    exit main
                end if
                x3 = x2
            endif
            l = 0
            if ( x1<=x2 ) then
                t2 = x1
                x1 = x2
                x2 = t2
                if ( t2<=x3 ) then
                    x2 = x3
                    x3 = t2
                endif
            endif
            x3 = x3 - t3
        end block main

        x1 = x1 - t3
        x2 = x2 - t3

        ! outputs:
        select case (l)
        case(0) ! three real roots
            zr = [x1,x2,x3]
            zi = 0.0_wp
        case(2) ! one real and two commplex roots
            zr = [x1, x2, x2]
            zi = [0.0_wp, x3, -x3]
        end select

    end subroutine rroots_chebyshev_cubic
!*****************************************************************************************

!*****************************************************************************************
!>
!  Sorts a set of complex numbers (with real and imag parts
!  in different vectors) in increasing order.
!
!  Uses a non-recursive quicksort, reverting to insertion sort on arrays of
!  size \(\le 20\). Dimension of `stack` limits array size to about \(2^{32}\).
!
!### License
!  * [Original LAPACK license](http://www.netlib.org/lapack/LICENSE.txt)
!
!### History
!  * Based on the LAPACK routine [DLASRT](http://www.netlib.org/lapack/explore-html/df/ddf/dlasrt_8f.html).
!  * Extensively modified by Jacob Williams,Feb. 2016. Converted to
!    modern Fortran and removed the descending sort option.

    pure subroutine sort_roots(x,y)

    implicit none

    real(wp),dimension(:),intent(inout) :: x  !! the real parts to be sorted.
                                              !! on exit,`x` has been sorted into
                                              !! increasing order (`x(1) <= ... <= x(n)`)
    real(wp),dimension(:),intent(inout) :: y  !! the imaginary parts to be sorted

    integer,parameter :: stack_size = 32 !! size for the stack arrays
    integer,parameter :: max_size_for_insertion_sort = 20 !! max size for using insertion sort.

    integer,dimension(2,stack_size) :: stack
    integer :: endd,i,j,n,start,stkpnt
    real(wp) :: d1,d2,d3
    real(wp) :: dmnmx,tmpx
    real(wp) :: dmnmy,tmpy

    ! number of elements to sort:
    n = size(x)

    if ( n>1 ) then

        stkpnt     = 1
        stack(1,1) = 1
        stack(2,1) = n

        do

            start  = stack(1,stkpnt)
            endd   = stack(2,stkpnt)
            stkpnt = stkpnt - 1
            if ( endd-start<=max_size_for_insertion_sort .and. endd>start ) then

                ! do insertion sort on x( start:endd )
                insertion: do i = start + 1,endd
                    do j = i,start + 1,-1
                        if ( x(j) < x(j-1)  ) then
                            dmnmx  = x(j)
                            x(j)   = x(j-1)
                            x(j-1) = dmnmx
                            dmnmy  = y(j)
                            y(j)   = y(j-1)
                            y(j-1) = dmnmy
                        else
                            exit
                        end if
                    end do
                end do insertion

            elseif ( endd-start>max_size_for_insertion_sort ) then

                ! partition x( start:endd ) and stack parts,largest one first
                ! choose partition entry as median of 3

                d1 = x(start)
                d2 = x(endd)
                i  = (start+endd)/2
                d3 = x(i)
                if ( d1 < d2 ) then
                    if ( d3 < d1 ) then
                        dmnmx = d1
                    elseif ( d3 < d2 ) then
                        dmnmx = d3
                    else
                        dmnmx = d2
                    endif
                elseif ( d3 < d2 ) then
                    dmnmx = d2
                elseif ( d3 < d1 ) then
                    dmnmx = d3
                else
                    dmnmx = d1
                endif

                i = start - 1
                j = endd + 1
                do
                    do
                        j = j - 1
                        if ( x(j) <= dmnmx ) exit
                    end do
                    do
                        i = i + 1
                        if ( x(i) >= dmnmx ) exit
                    end do
                    if ( i<j ) then
                        tmpx = x(i)
                        x(i) = x(j)
                        x(j) = tmpx
                        tmpy = y(i)
                        y(i) = y(j)
                        y(j) = tmpy
                    else
                        exit
                    endif
                end do
                if ( j-start>endd-j-1 ) then
                    stkpnt          = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                    stkpnt          = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                else
                    stkpnt          = stkpnt + 1
                    stack(1,stkpnt) = j + 1
                    stack(2,stkpnt) = endd
                    stkpnt          = stkpnt + 1
                    stack(1,stkpnt) = start
                    stack(2,stkpnt) = j
                endif

            endif

            if ( stkpnt<=0 ) exit

        end do

    end if

    ! check the imag parts:
    do i = 1, size(x)-1
        if (x(i)==x(i+1)) then
            if (y(i)>y(i+1)) then
                ! swap
                tmpy = y(i)
                y(i) = y(i+1)
                y(i+1) = tmpy
             end if
        end if
    end do

    end subroutine sort_roots
!*****************************************************************************************

!*****************************************************************************************
!>
!  Given coefficients `A(1),...,A(NDEG+1)` this subroutine computes the
!  `NDEG` roots of the polynomial `A(1)*X**NDEG + ... + A(NDEG+1)`
!  storing the roots as complex numbers in the array `Z`.
!  Require `NDEG >= 1` and `A(1) /= 0`.
!
!### Reference
!  * Original code from [JPL MATH77 Library](https://netlib.org/math/)
!
!### History
!  * C.L.Lawson & S.Y.Chan, JPL, June 3, 1986.
!  * 1987-09-16 DPOLZ  Lawson  Initial code.
!  * 1988-06-07 DPOLZ  CLL Reordered spec stmts for ANSI standard.
!  * 1988-11-16        CLL More editing of spec stmts.
!  * 1992-05-11 CLL IERR was not being set when N = 0 or 1. Fixed this. Added type stmts for all variables.
!  * 1992-05-11 DPOLZ  CLL
!  * 1994-10-19 DPOLZ  Krogh  Changes to use M77CON
!  * 1995-01-25 DPOLZ  Krogh Automate C conversion.
!  * 1995-11-17 DPOLZ  Krogh SFTRAN converted to Fortran 77
!  * 1996-03-30 DPOLZ  Krogh Added external statement, MIN0 => MIN.
!  * 1996-04-27 DPOLZ  Krogh Changes to use .C. and C%%.
!  * 2001-05-25 DPOLZ  Krogh Minor change for making .f90 version.
!  * 2022-09-24, Jacob Williams modernized this routine

    subroutine dpolz(ndeg,a,zr,zi,ierr)

        implicit none

        integer,intent(in) :: ndeg !! Degree of the polynomial
        real(wp),intent(in) :: a(ndeg+1) !! Contains the coefficients of a polynomial, high
                                         !! order coefficient first with `A(1)/=0`.
        real(wp),intent(out) :: zr(ndeg) !! Contains the real parts of the roots
        real(wp),intent(out) :: zi(ndeg) !! Contains the imaginary parts of the roots
        integer,intent(out) :: ierr !! Error flag:
                                    !!
                                    !! * Set by the subroutine to `0` on normal termination.
                                    !! * Set to `-1` if `A(1)=0`.
                                    !! * Set to `-2` if `NDEG<=0`.
                                    !! * Set to `J > 0` if the iteration count limit
                                    !!   has been exceeded and roots 1 through `J` have not been
                                    !!   determined.

        integer :: i,j,k,l,m,n,en,ll,mm,na,its,low,mp2,enm2
        real(wp) :: p,q,r,s,t,w,x,y,zz
        real(wp) :: c,f,g
        logical :: notlas , more
        real(wp),dimension(:,:),allocatable :: h !! Array of work space `(ndeg,ndeg)`
        real(wp),dimension(:),allocatable :: z !! Contains the polynomial roots stored as complex
                                               !! numbers. The real and imaginary parts of the Jth roots
                                               !! will be stored in `Z(2*J-1)` and `Z(2*J)` respectively.

        real(wp),parameter :: zero = 0.0_wp
        real(wp),parameter :: one  = 1.0_wp
        real(wp),parameter :: c75  = 0.75_wp
        real(wp),parameter :: half = 0.5_wp
        real(wp),parameter :: c43  = -0.4375_wp
        real(wp),parameter :: c95  = 0.95_wp
        real(wp),parameter :: machep = eps        !! d1mach(4)
        integer,parameter :: base = radix(1.0_wp) !! i1mach(10)
        integer,parameter :: b2 = base*base

        ierr = 0

        if ( ndeg<=0 ) then
            ierr = -2
            write(*,*) 'ndeg <= 0'
            return
        endif

        if ( a(1)==zero ) then
           ierr = -1
           write(*,*) 'a(1) == zero'
           return
        endif

        n = ndeg
        ierr = 0
        allocate(h(n,n)); h = zero ! workspace arrays
        allocate(z(2*n)); z = zero

        ! build first row of companion matrix.

        do i = 2 , ndeg + 1
           h(1,i-1) = -(a(i)/a(1))
        enddo

        ! extract any exact zero roots and set n = degree of
        ! remaining polynomial.

        do j = ndeg , 1 , -1
           if ( h(1,j)/=zero ) exit
           z(2*j-1) = zero
           z(2*j) = zero
           n = n - 1
        enddo

        ! special for n = 0 or 1.

        if ( n==0 ) return
        if ( n==1 ) then
           z(1) = h(1,1)
           return
        endif

        ! build rows 2 thru n of the companion matrix.

        do i = 2 , n
           do j = 1 , n
              h(i,j) = zero
           enddo
           h(i,i-1) = one
        enddo

  ! ***************** balance the matrix ***********************
  !
  !     this is an adaption of the eispack subroutine balanc to
  !     the special case of a companion matrix. the eispack balance
  !     is a translation of the algol procedure balance, num. math.
  !     13, 293-304(1969) by parlett and reinsch.
  !     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).

        do
            ! ********** iterative loop for norm reduction **********
            more = .false.
            do i = 1 , n
                ! compute r = sum of magnitudes in row i skipping diagonal.
                !         c = sum of magnitudes in col i skipping diagonal.
                if ( i==1 ) then
                    r = abs(h(1,2))
                    do j = 3 , n
                        r = r + abs(h(1,j))
                    enddo
                    c = abs(h(2,1))
                else
                    r = abs(h(i,i-1))
                    c = abs(h(1,i))
                    if ( i/=n ) c = c + abs(h(i+1,i))
                endif

                ! determine column scale factor, f.

                g = r/base
                f = one
                s = c + r

                do
                    if ( c>=g ) exit
                    f = f*base
                    c = c*b2
                end do
                g = r*base
                do
                    if ( c<g ) exit
                    f = f/base
                    c = c/b2
                end do

                ! will the factor f have a significant effect ?

                if ( (c+r)/f<c95*s ) then

                    ! yes, so do the scaling.

                    g = one/f
                    more = .true.

                    ! scale row i

                    if ( i==1 ) then
                        do j = 1 , n
                            h(1,j) = h(1,j)*g
                        enddo
                    else
                        h(i,i-1) = h(i,i-1)*g
                    endif

                    ! scale column i

                    h(1,i) = h(1,i)*f
                    if ( i/=n ) h(i+1,i) = h(i+1,i)*f

                endif
            enddo
            if ( .not. more ) exit
        end do

  ! ***************** qr eigenvalue algorithm ***********************
  !
  !     this is the eispack subroutine hqr that uses the qr
  !     algorithm to compute all eigenvalues of an upper
  !     hessenberg matrix. original algol code was due to martin,
  !     peters, and wilkinson, numer. math., 14, 219-231(1970).
  !
        low = 1
        en = n
        t = zero

        main : do

            ! ********** search for next eigenvalues **********
            if ( en<low ) exit main
            its = 0
            na = en - 1
            enm2 = na - 1

            sub : do
                ! ********** look for single small sub-diagonal element
                !            for l=en step -1 until low do -- **********
                do ll = low , en
                   l = en + low - ll
                   if ( l==low ) exit
                   if ( abs(h(l,l-1))<=machep*(abs(h(l-1,l-1))+abs(h(l,l))) ) exit
                enddo

                ! ********** form shift **********
                x = h(en,en)
                if ( l==en ) then
                   ! ********** one root found **********
                   z(2*en-1) = x + t
                   z(2*en) = zero
                   en = na
                else
                   y = h(na,na)
                   w = h(en,na)*h(na,en)
                   if ( l==na ) then
                      ! ********** two roots found **********
                      p = (y-x)*half
                      q = p*p + w
                      zz = sqrt(abs(q))
                      x = x + t
                      if ( q<zero ) then
                         ! ********** complex pair **********
                         z(2*na-1) = x + p
                         z(2*na) = zz
                         z(2*en-1) = x + p
                         z(2*en) = -zz
                      else
                         ! ********** pair of reals **********
                         zz = p + sign(zz,p)
                         z(2*na-1) = x + zz
                         z(2*na) = zero
                         z(2*en-1) = z(2*na-1)
                         z(2*en) = z(2*na)
                         if ( zz/=zero ) then
                            z(2*en-1) = x - w/zz
                            z(2*en) = zero
                         endif
                      endif
                      en = enm2
                   elseif ( its==30 ) then
                      ! ********** set error -- no convergence to an eigenvalue after 30 iterations **********
                      ierr = en
                      exit main
                   else
                      if ( its==10 .or. its==20 ) then
                         ! ********** form exceptional shift **********
                         t = t + x

                         do i = low , en
                            h(i,i) = h(i,i) - x
                         enddo

                         s = abs(h(en,na)) + abs(h(na,enm2))
                         x = c75*s
                         y = x
                         w = c43*s*s
                      endif
                      its = its + 1
                      !     ********** look for two consecutive small
                      !                sub-diagonal elements.
                      !                for m=en-2 step -1 until l do -- **********
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
                         if ( abs(h(m,m-1))*(abs(q)+abs(r))<=machep*abs(p) &
                              *(abs(h(m-1,m-1))+abs(zz)+abs(h(m+1,m+1))) ) &
                              exit
                      enddo

                      mp2 = m + 2

                      do i = mp2 , en
                         h(i,i-2) = zero
                         if ( i/=mp2 ) h(i,i-3) = zero
                      enddo
                      ! ********** double qr step involving rows l to en and
                      !            columns m to en **********
                      do k = m , na
                         notlas = k/=na
                         if ( k/=m ) then
                            p = h(k,k-1)
                            q = h(k+1,k-1)
                            r = zero
                            if ( notlas ) r = h(k+2,k-1)
                            x = abs(p) + abs(q) + abs(r)
                            if ( x==zero ) cycle !goto 640
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
                         ! ********** row modification **********
                         do j = k , en
                            p = h(k,j) + q*h(k+1,j)
                            if ( notlas ) then
                               p = p + r*h(k+2,j)
                               h(k+2,j) = h(k+2,j) - p*zz
                            endif
                            h(k+1,j) = h(k+1,j) - p*y
                            h(k,j) = h(k,j) - p*x
                         enddo

                         j = min(en,k+3)
                         ! ********** column modification **********
                         do i = l , j
                            p = x*h(i,k) + y*h(i,k+1)
                            if ( notlas ) then
                               p = p + zz*h(i,k+2)
                               h(i,k+2) = h(i,k+2) - p*r
                            endif
                            h(i,k+1) = h(i,k+1) - p*q
                            h(i,k) = h(i,k) - p
                         enddo

                      enddo
                      cycle sub
                   endif
                endif
                exit sub
            end do sub

        end do main

        if ( ierr/=0 ) write(*,*) 'convergence failure'

        ! return the computed roots:
        do i = 1, ndeg
            zr(i) = Z(2*i-1)
            zi(i) = Z(2*i)
        end do

    end subroutine dpolz
!*****************************************************************************************

!*****************************************************************************************
!>
!  In the discussion below, the notation A([*,],k} should be interpreted
!  as the complex number A(k) if A is declared complex, and should be
!  interpreted as the complex number A(1,k) + i * A(2,k) if A is not
!  declared to be of type complex.  Similar statements apply for Z(k).
!
!  Given complex coefficients A([*,[1),...,A([*,]NDEG+1) this
!  subr computes the NDEG roots of the polynomial
!              A([*,]1)*X**NDEG + ... + A([*,]NDEG+1)
!  storing the roots as complex numbers in the array Z( ).
!  Require NDEG >= 1 and A([*,]1) /= 0.
!
!### Reference
!  * Original code from [JPL MATH77 Library](https://netlib.org/math/)
!
!### License
!  Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
!  ALL RIGHTS RESERVED.
!  Based on Government Sponsored Research NAS7-03001.
!
!### History
!  * C.L.Lawson & S.Y.Chan, JPL, June 3,1986.
!  * 1987-02-25 CPOLZ  Lawson  Initial code.
!  * 1989-10-20 CLL Delcared all variables.
!  * 1992-05-11 CLL IERR was not being set when N = 0 or 1. Fixed this.
!  * 1995-01-18 CPOLZ  Krogh More M77CON for conversion to C.
!  * 1995-11-17 CPOLZ  Krogh Added M77CON statements for conversion to C
!  * 1995-11-17 CPOLZ  Krogh Converted SFTRAN to Fortran 77.
!  * 1996-03-30 CPOLZ  Krogh Added external statement.
!  * 1996-04-27 CPOLZ  Krogh Changes to use .C. and C%%.
!  * 2001-05-25 CPOLZ  Krogh Minor change for making .f90 version.
!  * 2022-10-06, Jacob Williams modernized this routine

    subroutine cpolz(a,ndeg,z,ierr)

    integer,intent(in) :: ndeg !! degree of the polynomial
    complex(wp),intent(in) :: a(ndeg+1) !! contains the complex coefficients of a polynomial
                                        !! high order coefficient first, with a([*,]1)/=0. the
                                        !! real and imaginary parts of the jth coefficient must
                                        !! be provided in a([*],j). the contents of this array will
                                        !! not be modified by the subroutine.
    complex(wp),intent(out) :: z(ndeg) !! contains the polynomial roots stored as complex
                                       !! numbers.  the real and imaginary parts of the jth root
                                       !! will be stored in z([*,]j).
    integer,intent(out) :: ierr !! error flag. set by the subroutine to 0 on normal
                                !! termination. set to -1 if a([*,]1)=0. set to -2 if ndeg
                                !! <= 0. set to  j > 0 if the iteration count limit
                                !! has been exceeded and roots 1 through j have not been
                                !! determined.

    complex(wp) :: temp
    integer :: i, j, n
    real(wp) :: c, f, g, r, s
    logical :: more, first
    real(wp) :: h(ndeg,ndeg,2) !! array of work space

    real(wp),parameter :: zero = 0.0_wp
    real(wp),parameter :: one = 1.0_wp
    real(wp),parameter :: c95 = 0.95_wp
    integer,parameter :: base = radix(1.0_wp)   !! i1mach(10)
    integer,parameter :: b2 = base * base

    if (ndeg <= 0) then
        ierr = -2
        write(*,*) 'ndeg <= 0'
        return
    end if

    if (a(1) == cmplx(zero, zero, wp)) then
        ierr = -1
        write(*,*) 'a(*,1) == zero'
        return
    end if

    n = ndeg
    ierr = 0

    ! build first row of companion matrix.

    do i = 2,n+1
        temp = -(a(i)/a(1))
        h(1,i-1,1) = real(temp,wp)
        h(1,i-1,2) = aimag(temp)
    end do

    ! extract any exact zero roots and set n = degree of
    ! remaining polynomial.

    do j = ndeg,1,-1
        if (h(1,j,1)/=zero .or. h(1,j,2)/=zero) exit
        z(j) = zero
        n = n - 1
    end do

    ! special for n = 0 or 1.

    if (n == 0) return
    if (n == 1) then
        z(1) = cmplx(h(1,1,1),h(1,1,2),wp)
        return
    end if

    ! build rows 2 thru n of the companion matrix.

    do i = 2,n
        do j = 1,n
            if (j == i-1) then
                h(i,j,1) = one
                h(i,j,2) = zero
            else
                h(i,j,1) = zero
                h(i,j,2) = zero
            end if
        end do
    end do

    ! ***************** balance the matrix ***********************
    !
    !     this is an adaption of the eispack subroutine balanc to
    !     the special case of a complex companion matrix. the eispack
    !     balance is a translation of the algol procedure balance,
    !     num. math. 13, 293-304(1969) by parlett and reinsch.
    !     handbook for auto. comp., vol.ii-linear algebra, 315-326(1971).

    ! ********** iterative loop for norm reduction **********
    do
        more = .false.
        do i = 1, n
          ! compute r = sum of magnitudes in row i skipping diagonal.
          !         c = sum of magnitudes in col i skipping diagonal.
          if (i == 1) then
            r = abs(h(1,2,1)) + abs(h(1,2,2))
            do j = 3,n
              r = r + abs(h(1,j,1)) + abs(h(1,j,2))
            end do
            c = abs(h(2,1,1)) + abs(h(2,1,2))
          else
            r = abs(h(i,i-1,1)) + abs(h(i,i-1,2))
            c = abs(h(1,i,1)) + abs(h(1,i,2))
            if (i /= n) then
              c = c + abs(h(i+1,i,1)) + abs(h(i+1,i,2))
            end if
          end if

          ! determine column scale factor, f.

          g = r / base
          f = one
          s = c + r

          do
            if (c >= g) exit
            f = f * base
            c = c * b2
          end do
          g = r * base
          do
            if (c < g) exit
            f = f / base
            c = c / b2
          end do
          ! will the factor f have a significant effect ?

          if ((c + r) / f < c95 * s) then

            ! yes, so do the scaling.

            g = one / f
            more = .true.

            ! scale row i

            if (i == 1) then
              do j = 1,n
                h(1,j,1) = h(1,j,1)*g
                h(1,j,2) = h(1,j,2)*g
              end do
            else
              h(i,i-1,1) = h(i,i-1,1)*g
              h(i,i-1,2) = h(i,i-1,2)*g
            end if

            ! scale column i

            h(1,i,1) = h(1,i,1) * f
            h(1,i,2) = h(1,i,2) * f
            if (i /= n) then
              h(i+1,i,1) = h(i+1,i,1) * f
              h(i+1,i,2) = h(i+1,i,2) * f
            end if

          end if
        end do
        if (.not. more) exit
    end do

    call scomqr(ndeg,n,1,n,h(1,1,1),h(1,1,2),z,ierr)

    if (ierr /= 0) write(*,*) 'Convergence failure in cpolz'

    end subroutine cpolz
!*****************************************************************************************

!*****************************************************************************************
!>
!  This subroutine finds the eigenvalues of a complex
!  upper hessenberg matrix by the qr method.
!
!  This subroutine is a translation of a unitary analogue of the
!  algol procedure  comlr, num. math. 12, 369-376(1968) by martin
!  and wilkinson.
!  handbook for auto. comp., vol.ii-linear algebra, 396-403(1971).
!  the unitary analogue substitutes the qr algorithm of francis
!  (comp. jour. 4, 332-345(1962)) for the lr algorithm.
!
!### Reference
!  * Original code from [JPL MATH77 Library](https://netlib.org/math/)
!
!### License
!  Copyright (c) 1996 California Institute of Technology, Pasadena, CA.
!  ALL RIGHTS RESERVED.
!  Based on Government Sponsored Research NAS7-03001.
!
!### History
!  * 1987-11-12 SCOMQR Lawson  Initial code.
!  * 1992-03-13 SCOMQR FTK  Removed implicit statements.
!  * 1995-01-03 SCOMQR WVS  Added EXTERNAL CQUO, CSQRT so VAX won't gripe
!  * 1996-01-18 SCOMQR Krogh  Added M77CON statements for conv. to C.
!  * 1996-03-30 SCOMQR Krogh  Added external statement.
!  * 1996-04-27 SCOMQR Krogh  Changes to use .C. and C%%.
!  * 2001-01-24 SCOMQR Krogh  CSQRT -> CSQRTX to avoid C lib. conflicts.
!  * 2022-10-06, Jacob Williams modernized this routine

    subroutine scomqr(nm,n,low,igh,hr,hi,z,ierr)

    integer,intent(in) :: nm !! the row dimension of two-dimensional array
                             !! parameters as declared in the calling program
                             !! dimension statement
    integer,intent(in) :: n !! the order of the matrix
    integer,intent(in) :: low !! low and igh are integers determined by the balancing
                              !! subroutine  cbal.  if  cbal  has not been used,
                              !! set low=1, igh=n
    integer,intent(in) :: igh !! low and igh are integers determined by the balancing
                              !! subroutine  cbal.  if  cbal  has not been used,
                              !! set low=1, igh=n
    real(wp),intent(inout) :: hi(nm,n) !! Input: hr and hi contain the real and imaginary parts,
                                       !! respectively, of the complex upper hessenberg matrix.
                                       !! their lower triangles below the subdiagonal contain
                                       !! information about the unitary transformations used in
                                       !! the reduction by  corth, if performed.
                                       !!
                                       !! Output: the upper hessenberg portions of hr and hi have been
                                       !! destroyed.  therefore, they must be saved before
                                       !! calling  comqr  if subsequent calculation of
                                       !! eigenvectors is to be performed,
    real(wp),intent(inout) :: hr(nm,n) !! see `hi` description
    complex(wp),intent(out) :: z(n) !! the real and imaginary parts,
                                    !! respectively, of the eigenvalues.  if an error
                                    !! exit is made, the eigenvalues should be correct
                                    !! for indices ierr+1,...,n,
    integer,intent(out) :: ierr !! is set to:
                                !!
                                !!  * zero -- for normal return,
                                !!  * j -- if the j-th eigenvalue has not been
                                !!    determined after 30 iterations.

    integer :: en,enm1,i,its,j,l,ll,lp1
    real(wp) :: norm,si,sr,ti,tr,xi,xr,yi,yr,zzi,zzr
    complex(wp) :: z3

    ierr = 0
    if (low /= igh) then
        ! create real subdiagonal elements
        l = low + 1

        do i = l, igh
            ll = min(i+1,igh)
            if (hi(i,i-1) == 0.0_wp) cycle
            norm = abs(cmplx(hr(i,i-1),hi(i,i-1),wp))
            yr = hr(i,i-1) / norm
            yi = hi(i,i-1) / norm
            hr(i,i-1) = norm
            hi(i,i-1) = 0.0_wp

            do j = i, igh
                si = yr * hi(i,j) - yi * hr(i,j)
                hr(i,j) = yr * hr(i,j) + yi * hi(i,j)
                hi(i,j) = si
            end do

            do j = low, ll
                si = yr * hi(j,i) + yi * hr(j,i)
                hr(j,i) = yr * hr(j,i) - yi * hi(j,i)
                hi(j,i) = si
            end do
        end do
    end if

      ! store roots isolated by cbal
      do i = 1, n
         if (i >= low .and. i <= igh) cycle
         z(i) = cmplx(hr(i,i),hi(i,i),wp)
      end do

      en = igh
      tr = 0.0_wp
      ti = 0.0_wp

    main : do
          ! search for next eigenvalue
          if (en < low) return
          its = 0
          enm1 = en - 1

          do
              ! look for single small sub-diagonal element
              ! for l=en step -1 until low
              do ll = low, en
                 l = en + low - ll
                 if (l == low) exit
                 if (abs(hr(l,l-1)) <= &
                          eps * (abs(hr(l-1,l-1)) + abs(hi(l-1,l-1)) &
                          + abs(hr(l,l)) +abs(hi(l,l)))) exit
              end do
              ! form shift
              if (l == en) then
              ! a root found
                z(en) = cmplx(hr(en,en)+tr,hi(en,en)+ti,wp)
                en = enm1
                cycle main
              end if
              if (its == 30) exit main
              if (its == 10 .or. its == 20) then
                ! form exceptional shift
                sr = abs(hr(en,enm1)) + abs(hr(enm1,en-2))
                si = 0.0_wp
              else
                sr = hr(en,en)
                si = hi(en,en)
                xr = hr(enm1,en) * hr(en,enm1)
                xi = hi(enm1,en) * hr(en,enm1)
                if (xr /= 0.0_wp .or. xi /= 0.0_wp) then
                    yr = (hr(enm1,enm1) - sr) / 2.0_wp
                    yi = (hi(enm1,enm1) - si) / 2.0_wp
                    z3 = sqrt(cmplx(yr**2-yi**2+xr,2.0_wp*yr*yi+xi,wp))
                    zzr = real(z3,wp)
                    zzi = aimag(z3)
                    if (yr * zzr + yi * zzi < 0.0_wp) then
                        zzr = -zzr
                        zzi = -zzi
                    end if
                    z3 = cmplx(xr,xi,wp) / cmplx(yr+zzr,yi+zzi,wp)
                    sr = sr - real(z3,wp)
                    si = si - aimag(z3)
                end if
              end if

              do i = low, en
                 hr(i,i) = hr(i,i) - sr
                 hi(i,i) = hi(i,i) - si
              end do

              tr = tr + sr
              ti = ti + si
              its = its + 1
              ! reduce to triangle (rows)
              lp1 = l + 1

              do i = lp1, en
                 sr = hr(i,i-1)
                 hr(i,i-1) = 0.0_wp
                 norm = sqrt(hr(i-1,i-1)*hr(i-1,i-1)+hi(i-1,i-1)*hi(i-1,i-1)+sr*sr)
                 xr = hr(i-1,i-1) / norm
                 xi = hi(i-1,i-1) / norm
                 z(i-1) = cmplx(xr,xi,wp)
                 hr(i-1,i-1) = norm
                 hi(i-1,i-1) = 0.0_wp
                 hi(i,i-1) = sr / norm
                 do j = i, en
                    yr = hr(i-1,j)
                    yi = hi(i-1,j)
                    zzr = hr(i,j)
                    zzi = hi(i,j)
                    hr(i-1,j) = xr * yr + xi * yi + hi(i,i-1) * zzr
                    hi(i-1,j) = xr * yi - xi * yr + hi(i,i-1) * zzi
                    hr(i,j) = xr * zzr - xi * zzi - hi(i,i-1) * yr
                    hi(i,j) = xr * zzi + xi * zzr - hi(i,i-1) * yi
                 end do
              end do

              si = hi(en,en)
              if (si /= 0.0_wp) then
                norm = abs(cmplx(hr(en,en),si,wp))
                sr = hr(en,en) / norm
                si = si / norm
                hr(en,en) = norm
                hi(en,en) = 0.0_wp
              end if
              ! inverse operation (columns)
              do j = lp1, en
                 xr = real(z(j-1),wp)
                 xi = aimag(z(j-1))
                 do i = l, j
                    yr = hr(i,j-1)
                    yi = 0.0
                    zzr = hr(i,j)
                    zzi = hi(i,j)
                    if (i /= j) then
                        yi = hi(i,j-1)
                        hi(i,j-1) = xr * yi + xi * yr + hi(j,j-1) * zzi
                    end if
                    hr(i,j-1) = xr * yr - xi * yi + hi(j,j-1) * zzr
                    hr(i,j) = xr * zzr + xi * zzi - hi(j,j-1) * yr
                    hi(i,j) = xr * zzi - xi * zzr - hi(j,j-1) * yi
                 end do
              end do

              if (si /= 0.0_wp) then
                do i = l, en
                    yr = hr(i,en)
                    yi = hi(i,en)
                    hr(i,en) = sr * yr - si * yi
                    hi(i,en) = sr * yi + si * yr
                end do
              end if

          end do

    end do main

    ! set error -- no convergence to an
    ! eigenvalue after 30 iterations
    ierr = en

    end subroutine scomqr

!*****************************************************************************************
end module polyroots_module
!*****************************************************************************************