submodule(stdlinalg) stdlinalg_qr
    implicit none

    contains
        ! Wrapper for calling LAPACK's real(dp) Q formation routine dorgqr.
        ! Forms the matrix Q in a QR or QRP factorization.
        ! The information about Q should be contained in Q, tau on input,
        ! and on output Q is changed to contain the full Q.
        module procedure qform_dp
            integer :: info

            call dorgqr(m, m, m, Q, m, tau, work, lwork, info)
        endprocedure qform_dp
endsubmodule stdlinalg_qr