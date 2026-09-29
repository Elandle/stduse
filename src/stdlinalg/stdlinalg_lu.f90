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


endsubmodule stdlinalg_lu