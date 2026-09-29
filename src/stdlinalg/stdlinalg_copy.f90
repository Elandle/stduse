submodule(stdlinalg) stdlinalg_copy
    implicit none

    contains
    
        !> \brief Sets \f$D = \text{diag}(A)\f$ for an \f$m \times n\f$ matrix \f$A\f$ and diagonal matrix \f$D\f$ stored as a vector.
        !!
        !! \param[inout] A (`real(dp), dimension(:, :)`) \f$n\times n\f$ matrix \f$A\f$ to update.
        !! \param[in]    D (`real(dp), dimension(min(size(A, 1), size(A, 2)))`)    Diagonal matrix \f$D\f$ stored as a vector.
        !! \param[in]    n (`integer`)                   Dimension of \f$A\f$ and \f$D\f$.
        module procedure diag_dp
            call diag_helper(A, D, size(A, 1), size(A, 2))
        contains
            subroutine diag_helper(A, D, m, n)
                integer , intent(in)  :: m
                integer , intent(in)  :: n
                real(dp), intent(in)  :: A(m, n)
                real(dp), intent(out) :: D(min(m, n))

                integer :: i

                do concurrent (i = 1 : min(m, n))
                    D(i) = A(i, i)
                enddo
            endsubroutine diag_helper
        endprocedure diag_dp

        !
        ! Sets:
        !
        ! B = uppertri(A)
        !
        ! That is, B is set to be the upper triangular part of A (including the diagonal),
        ! with all other entries set to zero.
        !
        module procedure uppertri_dp
            call uppertri_helper(A, B, size(A, 1), size(A, 2))
        contains
            subroutine uppertri_helper(A, B, m, n)
                integer , intent(in)  :: m
                integer , intent(in)  :: n
                real(dp), intent(in)  :: A(m, n)
                real(dp), intent(out) :: B(m, n)

                ! Zero the lower triangular part of B (including the diagonal)
                call dlaset('L', m, n, 0.0_dp, 0.0_dp, B, m)

                ! Copy the upper triangular part of A (including the diagonal) to B
                call dlacpy('U', m, n, A, m, B, m)
            endsubroutine uppertri_helper
        endprocedure uppertri_dp

        module procedure swap_dp
            real(dp) :: temp

            temp = x(i)
            x(i) = x(j)
            x(j) = temp
        endprocedure swap_dp

        ! Sets:
        !
        ! A = id
        !
        ! where id is the identity matrix and A is an n x n matrix
        module procedure make_identity_dp
            ! A = id
            call dlaset('A', n, n, 0.0_dp, 1.0_dp, A, n)
        endprocedure make_identity_dp

        ! Sets:
        !
        ! A = 0
        module procedure zero_matrix_dp
            ! A = 0
            call dlaset('A', n, n, 0.0_dp, 0.0_dp, A, n)
        endprocedure zero_matrix_dp

        !> Copies \f$A = B\f$, where \f$A\f$ and \f$B\f$ are both \f$m\times n\f$ matrices.
        !!
        !! \param[out]  A  (`real(dp), dimension(:, :)`) Matrix to copy into.
        !! \param[in]   B  (`real(dp), dimension(size(A, 1), size(A, 2))`) Matrix to copy from.
        module procedure copy_matrix_dp
            call copy_matrix_helper(A, B, size(A, 1), size(A, 2))
            contains
                subroutine copy_matrix_helper(A, B, m, n)
                    integer , intent(in)  :: m
                    integer , intent(in)  :: n
                    real(dp), intent(out) :: A(m, n)
                    real(dp), intent(in)  :: B(m, n)

                    ! A = B
                    call dlacpy('A', m, n, B, m, A, m)
                endsubroutine copy_matrix_helper
        endprocedure copy_matrix_dp

        ! Sets A to the transpose of B.
        module procedure transpose_dp
            integer :: i, j, m, n

            m = size(A, 1) ; n = size(A, 2)

            ! I believe Intel's MKL has something for this.

            do j = 1, n
                do i = 1, n
                    A(i, j) = B(j, i)
                enddo
            enddo
        endprocedure transpose_dp

        module procedure extractdiag_dp
            integer :: i

            do i = 1, m
                D(i) = A(i, i)
            enddo
        endprocedure extractdiag_dp

        module procedure id_dp
            call dlaset('A', m, m, 0.0_dp, 1.0_dp, A, m)
        endprocedure id_dp






endsubmodule stdlinalg_copy