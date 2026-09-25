submodule(stdlinalg) stdlinalg_permutation
    implicit none

    contains

        module procedure invert_permutation
            ! Knuth, The Art of Computer Programming, Volume 1
            ! Section 1.3.3, algorithm I
            integer :: m, n, j, i

            n = size(P) ; if (n .eq. 0) return

            ! I1
            m = n
            j = -1

            ! I2
20          i = P(m)
            if (i .lt. 0) then
                goto 50
            endif

            ! I3
30          P(m) = j
            j = -m
            m = i
            i = P(m)

            ! I4
            if (i .gt. 0) then
                goto 30
            else
                i = j
            endif

            ! I5
50          P(m) = -i

            ! I6
            m = m - 1
            if (m .gt. 0) then
                goto 20
            endif
        endprocedure invert_permutation

        module procedure colpivswap_dp
            !
            ! Swaps the columns of the n x n matrix A according to the n
            ! long integer vector piv
            !
            ! If piv(i) = k, then column i of A becomes column k of A
            !
            ! In other words, for i = 1, 2, ..., n:
            !
            !       A(:, piv(i)) = A(:, i)
            !
            ! Where this replacement is done independently (changes where
            ! columns are do not affect other column changes)
            !

            integer i

            !
            ! TODO:
            ! Can this be done without copying A?
            !

            call copy_matrix(matwork, A)

            do i = 1, n
               if (piv(i) .ne. i) then ! Do not copy a column if it is already in the right place
                   ! A(:, piv(i)) = matwork(:, i)
                   call dcopy(n, matwork(1, i), 1, A(1, piv(i)), 1)
               endif
            enddo
        endprocedure colpivswap_dp








endsubmodule stdlinalg_permutation