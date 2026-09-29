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

        ! Sets:
        !
        ! I = inv(P)
        !
        ! In terms of permutation matrices stored as length n vectors.
        ! Todo: this indeed inverts a permutation P, but not the ipiv type
        ! in LAPACK and doesn't do any checking. Clarify this usage and remove the old tag.
        module procedure invert_permutation_old
            integer :: j

            I(P) = [(j, j = 1, size(P))]
        endprocedure invert_permutation_old

        ! Updates:
        !
        ! A = A * P
        !
        ! Where P is a permutation matrix stored as a vector.
        !
        ! This subroutine is implemented basically as a copy of Julia's intrinsic permutecols
        module procedure permutecols_dp
            integer :: count
            integer :: start
            integer :: ptr
            integer :: next

            count = 0; start = 0
        
            do while (count .lt. n)
                start = findloc(P, start + 1, mask = P .ne. 0, dim = 1)
                if (start .eq. 0) exit
                ptr = start

                ! Equivalent (without using findloc):
                ! do start = start + 1, n
                !     if (P(start) .ne. 0) then
                !         ptr = start
                !         exit
                !     endif
                ! enddo
        
                next = P(start)
                count = count + 1
        
                do while (next .ne. start)
                    call dswap(n, A(1, ptr), 1, A(1, next), 1) ! Swap columns A(:, ptr) <--> A(:, next)
                    P(ptr) = 0
                    ptr = next
                    next = P(next)
                    count = count + 1
                end do

                P(ptr) = 0
            end do
        endprocedure permutecols_dp

        ! These following permutation algorithms are all based off the same underlying algorithm
        ! TODO:
        ! Find an inplace permutation algorithm that doesn't need an inverse permutation
        module procedure permute_matrix_columns_dp
            integer :: i, j, k

            do i = 1, m
                j = P(i)
                if (j .ne. i) then
                    k = Pinv(i)
                    do
                        if (i .lt. j .and. i .lt. k) then
                            if (j .eq. k) then
                                call rotate_matrix_columns(P, i, A, m)
                                exit
                            endif
                            j = P(j)
                            if (j .eq. k) then
                                call rotate_matrix_columns(P, i, A, m)
                                exit
                            endif
                            k = Pinv(k)
                        else
                            exit
                        endif
                    enddo
                endif
            enddo
        endprocedure permute_matrix_columns_dp

        module procedure permute_matrix_rows_dp
            integer :: i, j, k

            do i = 1, m
                j = P(i)
                if (j .ne. i) then
                    k = Pinv(i)
                    do
                        if (i .lt. j .and. i .lt. k) then
                            if (j .eq. k) then
                                call rotate_matrix_rows(P, i, A, m)
                                exit
                            endif
                            j = P(j)
                            if (j .eq. k) then
                                call rotate_matrix_rows(P, i, A, m)
                                exit
                            endif
                            k = Pinv(k)
                        else
                            exit
                        endif
                    enddo
                endif
            enddo
        endprocedure permute_matrix_rows_dp

        ! Part of permute
        module procedure rotate_matrix_columns_dp
            integer :: i

            i = P(leader)
            do
                if (i .eq. leader) then
                    exit
                else
                    call dswap(m, A(1, i), 1, A(1, leader), 1)
                    i = P(i)
                endif
            enddo
        endprocedure rotate_matrix_columns_dp

        ! Part of permute
        module procedure rotate_matrix_rows_dp
            integer :: i

            i = P(leader)
            do
                if (i .eq. leader) then
                    exit
                else
                    call dswap(m, A(i, 1), m, A(leader, 1), m)
                    i = P(i)
                endif
            enddo
        endprocedure rotate_matrix_rows_dp

        ! Fich, Munro, and Poblete
        ! Permuting in Place
        ! Figure 5
        module procedure permute_dp
            integer :: i, j, k

            do i = 1, m
                j = P(i)
                if (j .ne. i) then
                    k = Pinv(i)
                    do
                        if (i .lt. j .and. i .lt. k) then
                            if (j .eq. k) then
                                call rotate(P, i, x)
                                exit
                            endif
                            j = P(j)
                            if (j .eq. k) then
                                call rotate(P, i, x)
                                exit
                            endif
                            k = Pinv(k)
                        else
                            exit
                        endif
                    enddo
                endif
            enddo
        endprocedure permute_dp

        ! Part of permute
        module procedure rotate_dp
            integer :: i

            i = P(leader)
            do
                if (i .eq. leader) then
                    exit
                else
                    call swap(x, i, leader)
                    i = P(i)
                endif
            enddo
        endprocedure rotate_dp

        module procedure dlaswpc
            !    MODIFIED VERSION OF DLASWP
            !  -- LAPACK auxiliary routine --
            !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
            !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
            ! =====================================================================
            !
            !     .. Local Scalars ..
            integer            i, i1, i2, inc, ip, ix, ix0, j, k, n32
            double precision   temp
            !     ..
            !     .. Executable Statements ..
            !
            !     Interchange row i with row ipiv(k1+(i-k1)*abs(incx)) for each of rows
            !     k1 through k2.
            !
            if (incx .gt. 0) then
                ix0 = k1
                i1  = k1
                i2  = k2
                inc = 1
            elseif (incx .lt. 0) then
                ix0 = k1 + (k1-k2) * incx
                i1  = k2
                i2  = k1
                inc = -1
            else
                return
            endif
            !
            n32 = (n/32) * 32
            if (n32 .ne. 0) then
                do 30 j = 1, n32, 32
                    ix = ix0
                    do 20 i = i1, i2, inc
                        ip = ipiv(ix)
                        if (ip .ne. i) then
                            do 10 k = j, j + 31
                                temp     = a(k , i)
                                a(k , i) = a(k, ip)
                                a(k, ip) = temp
        10                  continue
                        endif
                        ix = ix + incx
        20          continue
        30      continue
            endif
        
            if (n32 .ne. n) then
                n32 = n32 + 1
                ix  = ix0
                do 50 i = i1, i2, inc
                    ip = ipiv(ix)
                    if (ip .ne. i) then
                        do 40 k = n32, n
                            temp     = a(k , i)
                            a(k , i) = a(k, ip)
                            a(k, ip) = temp
        40              continue
                    endif
                    ix = ix + incx
        50      continue
            endif
            !
            return
            !
            !     End of DLASWP
            !
        endprocedure dlaswpc

        ! Computes the sign of a piv permutation.
        module procedure pivsgn
            integer :: apivl, i

            if (present(pivl)) then
                apivl = pivl
            else
                apivl = size(piv)
            endif

            if ((apivl .gt. size(piv)) .or. (apivl .lt. 0)) then
                error stop "invalid pivl"
            endif

            sgn = 1
            do i = 1, apivl
                if (piv(i) .ne. i) sgn = -sgn
            enddo
        endprocedure pivsgn








endsubmodule stdlinalg_permutation