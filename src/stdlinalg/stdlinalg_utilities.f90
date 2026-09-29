submodule(stdlinalg) stdlinalg_utilities
    implicit none
    contains
        ! Returns the 2 norm of the matrix A
        module procedure twonorm_dp
            real(dp), allocatable :: B(:, :)
            real(dp), allocatable :: S(:)
            real(dp), allocatable :: work(:)
            integer               :: m
            integer               :: n
            integer               :: lwork
            integer               :: info
            
            m = size(A, 1) ; n = size(A, 2)
            lwork = 5 * (3 * min(m, n) + max(m, n) + 5 * min(m, n))

            allocate(B(m, n))
            allocate(S(min(m, n)))
            allocate(work(lwork))

            call dlacpy('A', m, n, A, m, B, m)
            call dgesvd('N', 'N', m, n, B, m, S, work, lwork, work, lwork, work, lwork, info)

            norm = S(1) 
        endprocedure twonorm_dp

        ! Updates:
        !
        ! A = A + transpose(B)
        !
        ! where A is m x n and B is n x m.
        module procedure add_transpose_dp
            call add_transpose_helper(A, B, size(A, 1), size(A, 2))
            contains
                subroutine add_transpose_helper(A, B, m, n)
                    integer , intent(in)    :: m
                    integer , intent(in)    :: n
                    real(dp), intent(inout) :: A(m, n)
                    real(dp), intent(in)    :: B(n, m)

                    integer :: j

                    ! Intel's MKL might have something that does this, but this is more portable.
                    do j = 1, n
                        call daxpy(m, 1.0_dp, B(j, 1), size(B, 1), A(1, j), 1)
                    enddo
                endsubroutine add_transpose_helper
        endprocedure add_transpose_dp

        ! Updates:
        !
        ! A = A + B
        !
        ! where A and B are m x n matrices
        module procedure add_matrix_dp
            call add_matrix_helper(A, B, size(A, 1), size(A, 2))
            contains
                subroutine add_matrix_helper(A, B, m, n)
                    integer , intent(in)    :: m
                    integer , intent(in)    :: n
                    real(dp), intent(inout) :: A(m, n)
                    real(dp), intent(in)    :: B(m, n)

                    call daxpy(size(A, 1) * size(A, 2), 1.0_dp, B, 1, A, 1)
                endsubroutine add_matrix_helper
        endprocedure add_matrix_dp

        module procedure avgdiag_dp
            integer  :: j

            avg = 0.0_dp
            do j = 1, m
                avg = avg + A(j, j)
            enddo
            avg = avg / m
        endprocedure avgdiag_dp

        module procedure matdiff_dp
            character(:), allocatable :: actualdifftype

            if (.not. present(difftype)) then
                actualdifftype = "difflim"
            else
                actualdifftype = difftype
            endif

            work = A - B
            ! Second work is never referenced in calling dlange (it is
            ! only needed when norm = 'I')
            if     (actualdifftype .eq. "difflim") then
                diff = dlange('F', m, n, work, m, work) / (m * n)
            elseif (actualdifftype .eq. "maxabs") then
                diff = dlange('M', m, n, work, m, work)
            elseif (actualdifftype .eq. "frobenius") then
                diff = dlange('F', m, n, work, m, work)
            elseif (actualdifftype .eq. "abssum") then
                diff = sum(abs(work))
            endif
        endprocedure matdiff_dp

        ! Effectively interchanges the names of the allocatable real(dp) matrices A and B.
        ! If A is matrix m1 and B matrix m2, then after calling this
        ! A is matrix m2 and B is matrix m1.
        module procedure alloc_matrixswap_dp
            real(dp), allocatable :: temp(:, :)

            call move_alloc(A, temp)
            call move_alloc(B, A)
            call move_alloc(temp, B)
        endprocedure alloc_matrixswap_dp

        !> Appends a column to the @c real(dp) matrix @p A
        !!
        !! @param[in,out] A Matrix to append a column to.
        !! @param[in]     x Vector to append to @p A.
        subroutine append_column_dp(A, x)
            real(dp), allocatable, intent(inout) :: A(:, :)
            real(dp),              intent(in)    :: x(:)

            real(dp), allocatable :: temp(:, :)

            ! Make x the first column of A if A is not already allocated
            if (.not. allocated(A)) then
                allocate(A(size(x), 1))
                A(:, 1) = x
                return
            endif

            ! Make sure A and x agree in dimension.
            if (size(A, 1) .ne. size(x)) error stop "error stop in procedure append_column_dp from module lattice_mod: attempting to append a column to a matrix with mismatching number of rows."

            allocate(temp(size(x), size(A, 2) + 1))
            temp(:, 1:size(A, 2)  ) = A
            temp(:, size(A, 2) + 1) = x
            call move_alloc(temp, A)
        endsubroutine append_column_dp



endsubmodule stdlinalg_utilities