submodule(stdlinalg) stdlinalg_multiply
    implicit none
    contains
        !> \brief Updates \f$A = AD\f$ for an \f$m \times n\f$ matrix \f$A\f$ and \f$m \times m\f$ diagonal matrix \f$D\f$ stored as a vector.
        !!
        !! \param[inout] A (`real(dp), dimension(:, :)`) \f$m\times n\f$ matrix \f$A\f$ to update.
        !! \param[in]    D (`real(dp), dimension(size(A, 1))`)    Diagonal matrix \f$D\f$ stored as a vector.
        module procedure right_diagmult_dp
            call right_diagmult_helper(A, D, size(A, 1), size(A, 2))
        contains
            subroutine right_diagmult_helper(A, D, m, n)
                integer , intent(in)    :: m
                integer , intent(in)    :: n
                real(dp), intent(inout) :: A(m, n)
                real(dp), intent(in)    :: D(n)

                integer :: i

                ! Scale column i of A by D(i)
                do i = 1, n
                    call dscal(m, D(i), A(1, i), 1)
                enddo
            endsubroutine right_diagmult_helper
        endprocedure right_diagmult_dp

        !> \brief Updates \f$A = DA\f$ for an \f$m \times n\f$ matrix \f$A\f$ and \f$n \times n\f$ diagonal matrix \f$D\f$ stored as a vector.
        !!
        !! \param[inout] A (`real(dp), dimension(:, :)`) \f$m\times n\f$ matrix \f$A\f$ to update.
        !! \param[in]    D (`real(dp), dimension(size(A, 2))`)    Diagonal matrix \f$D\f$ stored as a vector.
        !! \see Some BLAS/LAPACK distributions contain the `dlascl2` subroutine and some do not.
        !! `dlascl2` should be used if possible. If not, provided alternative code can be used.
        module procedure left_diagmult_dp
            call left_diagmult_helper(A, D, size(A, 1), size(A, 2))
        contains
            subroutine left_diagmult_helper(A, D, m, n)
                integer , intent(in)    :: m
                integer , intent(in)    :: n
                real(dp), intent(inout) :: A(m, n)
                real(dp), intent(in)    :: D(m)

                ! I have ran across BLAS/LAPACK's without this routine for some reason.
                ! If this routine is missing, comment out this line and uncomment the following four line
                ! integer declaration and for loop.
                call dlascl2(m, n, D, A, m)

                ! Scale row i of A by D(i)
                ! UNCOMMENT THESE FOUR LINES ------------
                ! integer :: i
                ! do i = 1, m
                !     call dscal(n, D(i), A(i, 1), m)
                ! enddo
                ! ----------------------------------------
            endsubroutine left_diagmult_helper
        endprocedure left_diagmult_dp

        !> \brief Updates \f$A = AD^{-1}\f$ for an \f$m \times n\f$ matrix \f$A\f$ and \f$n \times n\f$ diagonal matrix \f$D\f$ stored as a vector.
        !!
        !! \param[inout] A (`real(dp), dimension(:, :)`) \f$m\times n\f$ matrix \f$A\f$ to update.
        !! \param[in]    D (`real(dp), dimension(size(A, 2))`)    Diagonal matrix \f$D\f$ stored as a vector.
        module procedure right_diaginvmult_dp
            call right_diaginvmult_helper(A, D, size(A, 1), size(A, 2))
        contains
            subroutine right_diaginvmult_helper(A, D, m, n)
                integer , intent(in)    :: m
                integer , intent(in)    :: n
                real(dp), intent(inout) :: A(m, n)
                real(dp), intent(in)    :: D(n)

                integer :: i

                ! Scale column i of A by 1/D(i)
                do i = 1, n
                    call dscal(m, 1.0_dp / D(i), A(1, i), 1)
                enddo
            endsubroutine right_diaginvmult_helper
        endprocedure right_diaginvmult_dp

        !> \brief Updates \f$A = D^{-1}A\f$ for an \f$m \times n\f$ matrix \f$A\f$ and \f$m \times m\f$ diagonal matrix \f$D\f$ stored as a vector.
        !!
        !! \param[inout] A (`real(dp), dimension(:, :)`) \f$m\times n\f$ matrix \f$A\f$ to update.
        !! \param[in]    D (`real(dp), dimension(size(A, 2))`)    Diagonal matrix \f$D\f$ stored as a vector.
        module procedure left_diaginvmult_dp
            call left_diaginvmult_helper(A, D, size(A, 1), size(A, 2))
        contains
            subroutine left_diaginvmult_helper(A, D, m, n)
                integer , intent(in)    :: m
                integer , intent(in)    :: n
                real(dp), intent(inout) :: A(m, n)
                real(dp), intent(in)    :: D(m)

                integer :: i

                ! Scale row i of A by 1/D(i)
                do i = 1, m
                    call dscal(n, 1.0_dp / D(i), A(i, 1), m)
                enddo
            endsubroutine left_diaginvmult_helper
        endprocedure left_diaginvmult_dp

        ! Updates:
        !
        ! A = B * A
        !
        ! where A is an m x n matrix and B an m x m matrix (so the result can still be stored in A).
        ! work is a workspace array (any dimensional, eg 1 or 2 as long as it has enough space)
        ! at least m * n long.
        module procedure left_matmul_dp
            call left_matmul_helper(A, B, work, size(A, 1), size(A, 2))
        contains
            subroutine left_matmul_helper(A, B, work, m, n)
                integer , intent(in)    :: m
                integer , intent(in)    :: n
                real(dp), intent(inout) :: A(m, n)
                real(dp), intent(in)    :: B(m, m)
                real(dp), intent(out)   :: work(m*n)

                ! work = A
                call dlacpy('a', m, n, A, m, work, m)

                ! A = B * work
                call dgemm('n', 'n', m, n, m, 1.0_dp, B, m, work, m, 0.0_dp, A, m)
            endsubroutine left_matmul_helper
        endprocedure left_matmul_dp

        ! Updates:
        !
        ! A = A * B
        !
        ! where A is an m x n matrix and B an n x n matrix (so the result can still be stored in A).
        ! work is a workspace array (any dimensional, eg 1 or 2 as long as it has enough space)
        ! at least m * n long.
        module procedure right_matmul_dp
            call right_matmul_helper(A, B, work, size(A, 1), size(A, 2))
        contains
            subroutine right_matmul_helper(A, B, work, m, n)
                integer , intent(in)    :: m
                integer , intent(in)    :: n
                real(dp), intent(inout) :: A(m, n)
                real(dp), intent(in)    :: B(n, n)
                real(dp), intent(out)   :: work(m*n)

                ! work = A
                call dlacpy('a', m, n, A, m, work, m)

                ! A = work * B
                call dgemm('n', 'n', m, n, n, 1.0_dp, work, m, B, n, 0.0_dp, A, m)
            endsubroutine right_matmul_helper
        endprocedure right_matmul_dp


endsubmodule stdlinalg_multiply