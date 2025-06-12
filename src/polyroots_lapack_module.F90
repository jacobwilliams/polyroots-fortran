!*****************************************************************************************
!>
!  LAPACK interface module
!
!### Author
!  * Jacob Williams

module polyroots_lapack_module

#ifdef USE_STDLIB_LAPACK
#ifdef REAL32
    use stdlib_linalg_lapack, sgeev => stdlib_sgeev, cgeev => stdlib_cgeev
#elif REAL128
    use stdlib_linalg_lapack, qgeev => stdlib_qgeev, wgeev => stdlib_wgeev
#else
    use stdlib_linalg_lapack, dgeev => stdlib_dgeev, zgeev => stdlib_zgeev
#endif
#endif

implicit none

#ifndef USE_STDLIB_LAPACK

#ifdef REAL32
    interface
        subroutine sgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
            implicit none
            character :: jobvl, jobvr
            integer :: info, lda, ldvl, ldvr, lwork, n
            real :: a(lda, *), vl(ldvl, *), vr(ldvr, *), wi(*), work(*), wr(*)
        end subroutine sgeev
        subroutine cgeev( jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info )
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
        subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
            implicit none
            character :: jobvl, jobvr
            integer :: info, lda, ldvl, ldvr, lwork, n
            double precision :: a(lda, *), vl(ldvl, *), vr(ldvr, *), wi(*), work(*), wr(*)
        end subroutine dgeev
        subroutine zgeev( jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info )
            implicit none
            character :: jobvl, jobvr
            integer :: info, lda, ldvl, ldvr, lwork, n
            double precision :: rwork( * )
            complex*16 :: a( lda, * ), vl( ldvl, * ), vr( ldvr, * ), w( * ), work( * )
        end subroutine zgeev
    end interface
#endif
#endif

end module polyroots_lapack_module