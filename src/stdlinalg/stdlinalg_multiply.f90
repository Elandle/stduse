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


endsubmodule stdlinalg_multiply