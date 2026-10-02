submodule(stdlinalg) stdlinalg_lu
    implicit none
    contains
        ! Sets:
        !
        ! A = inv(A)
        !
        ! where A is a general matrix.
        !
        ! A is inverted in two steps:
        !
        ! First, A is A = P * L * U factorised by dgetrf
        ! Then, A is inverted by dgetri (using the result of dgetrf)
        !
        module procedure invert_dp
            integer :: i

            ! A = P * L * U
            call dgetrf(n, n, A, n, P, info)

            ! Calculation of the sign of the determinant
            if (present(detsgn)) then
                detsgn = 1
                do i = 1, n
                    if (A(i, i) .lt. 0.0_dp) then
                        detsgn = -detsgn
                    endif
                    if (P(i) .ne. i) then
                        detsgn = -detsgn
                    endif
                enddo
            endif

            ! A = inv(A)
            call dgetri(n, A, n, P, work, lwork, info)
        endprocedure invert_dp

        ! Sets:
        !
        ! det = det(A)
        !
        ! by first factorising A = P * L * U (destroying A),
        ! and calculating the determinant by multiplying the diagonal of U
        ! and adjusting the sign based on P
        module procedure determinant_dp

            integer :: i

            ! A = P * L * U
            call dgetrf(n, n, A, n, P, info)

            det = A(1, 1)
            if (P(1) .ne. 1) then
                det = -1 * det
            endif

            do i = 2, n
                det = det * A(i, i)
                if (P(i) .ne. i) then
                    det = -1 * det
                endif
            enddo
        endprocedure determinant_dp

        ! Sets:
        !
        ! det = det(A)
        !
        ! by first factorising A = P * L * U (destroying A),
        ! and calculating the determinant by multiplying the diagonal of U
        ! and adjusting the sign based on P
        module procedure determinant_zp
            integer :: i

            ! A = P * L * U
            call zgetrf(n, n, A, n, P, info)

            det = A(1, 1)
            if (P(1) .ne. 1) then
                det = -1 * det
            endif

            do i = 2, n
                det = det * A(i, i)
                if (P(i) .ne. i) then
                    det = -1 * det
                endif
            enddo
        endprocedure determinant_zp

        module procedure lu_dp
            integer :: m, n, lpiv, info

            m = size(A, 1) ; n = size(A, 2) ; lpiv = size(piv) ; info = 0

            ! Make sure piv can hold all pivots.
            if (lpiv .lt. min(m, n)) then
                error stop "error stop in procedure lu_dp from submodule stdlinalg_lu of module &
                           &linalg: piv array must have length at least min(size(A, 1), size(A, 2))."
            endif

            ! Make sure A is actually a matrix.
            if ((m .eq. 0) .or. (n .eq. 0)) return

            ! PA = LU factorize A.
            call dgetrf(m, n, A, m, piv, info)

            ! Check for an LAPACK error.
            if (info .ne. 0) then
                block
                    character(len=128) :: errmsg
                    write(errmsg, "(a, i0)") "lapack dgetrf call failed with info = ", info
                    error stop trim(errmsg)
                endblock
            endif
        endprocedure lu_dp

        module procedure ludet_dp
            integer  :: m, n, i

            m = size(A, 1) ; n = size(A, 2)

            if (m .ne. n) then
                error stop "input is not square."
            endif

            if (size(piv) < m) then
                error stop "piv array must have length at least the same size as A."
            end if

            det = real(pivsgn(piv, n), dp)
            do i = 1, n
                det = det * A(i, i)
            enddo
        endprocedure ludet_dp

        module procedure luinv_dp
            integer :: m, n, lwork, info

            m = size(A, 1) ; n = size(A, 2) ; lwork = size(work) ; info = 0

            ! Make sure A is square.
            if (m .ne. n) then
                error stop "input matrix is not square."
            end if

            ! Make sure piv length matches the dimension of A.
            if (size(piv) .lt. m) then
                error stop "piv dimension is not at least as long as the dimension of A."
            end if

            ! Make sure matrix is not trivial.
            if (m .eq.0) return

            ! dgetri requires a workspace at least m long.
            if (lwork .lt. m) then
                error stop "work array is not long enough."
            end if

            ! Invert the PA = LU factorization stored in A and piv.
            call dgetri(m, A, m, piv, work, lwork, info)

            ! Check for an LAPACK error.
            if (info .ne. 0) then
                block
                    character(len=128) :: errmsg
                    write(errmsg, "(a, i0)") "lapack dgetri call failed with info = ", info
                    error stop trim(errmsg)
                endblock
            endif
        end procedure luinv_dp


endsubmodule stdlinalg_lu