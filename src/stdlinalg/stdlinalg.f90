module stdlinalg
    use stduse
    implicit none

    interface right_diagmult
        module subroutine right_diagmult_dp(A, D)
            real(dp), contiguous, intent(inout) :: A(:, :)
            real(dp), contiguous, intent(in)    :: D(:)
        endsubroutine right_diagmult_dp
    endinterface right_diagmult

    interface left_diagmult
        module subroutine left_diagmult_dp(A, D)
            real(dp), contiguous, intent(inout) :: A(:, :)
            real(dp), contiguous, intent(in)    :: D(:)
        endsubroutine left_diagmult_dp
    endinterface left_diagmult

    interface right_diaginvmult
        module subroutine right_diaginvmult_dp(A, D)
            real(dp), contiguous, intent(inout) :: A(:, :)
            real(dp), contiguous, intent(in)    :: D(:)
        endsubroutine right_diaginvmult_dp
    endinterface right_diaginvmult

    interface left_diaginvmult
        module subroutine left_diaginvmult_dp(A, D)
            real(dp), contiguous, intent(inout) :: A(:, :)
            real(dp), contiguous, intent(in)    :: D(:)
        endsubroutine left_diaginvmult_dp
    endinterface left_diaginvmult


    interface twonorm
        module function twonorm_dp(A) result(norm)
            real(dp), intent(in)  :: A(:, :)
            real(dp) :: norm
        endfunction twonorm_dp
    endinterface twonorm


    ! stdlinalg_matexp submodule interfaces ------------------------------
    interface gpadm
        module subroutine dgpadm(ideg, m, t, H, ldh, wsp, lwsp, ipiv, iexph, ns, iflag)
            integer  :: ideg, m, ldh, lwsp, iexph, ns, iflag, ipiv(m)
            real(dp) :: t, H(ldh,m), wsp(lwsp)
        endsubroutine dgpadm
    endinterface gpadm
    
    interface expm
        module subroutine dexpm(A, exptA, t, ideg)
            real(dp), intent(in)  :: A(:, :)
            real(dp), intent(out) :: exptA(size(A, 1), size(A, 1))
            real(dp), optional    :: t
            integer , optional    :: ideg
        endsubroutine dexpm
    endinterface expm
    ! ---------------------------------------------------------------------

    ! stdlinalg_permutations submodule interfaces ------------------------

    interface invert_permutation
        module subroutine invert_permutation(P)
            integer, intent(inout) :: P(:)
        endsubroutine invert_permutation
    endinterface

    interface colpivswap
        module subroutine colpivswap_dp(A, piv, n, matwork)
            real(dp), intent(inout) :: A(n, n)
            integer , intent(in)    :: piv(n)
            integer , intent(in)    :: n
            real(dp), intent(out)   :: matwork(n, n)
        endsubroutine colpivswap_dp
    endinterface

    ! ---------------------------------------------------------------------

        interface diag
            module subroutine diag_dp(A, D)
                real(dp), contiguous, intent(in)  :: A(:, :)
                real(dp), contiguous, intent(out) :: D(:)
            endsubroutine diag_dp
        endinterface diag

        interface uppertri
            module subroutine uppertri_dp(A, B)
                real(dp), contiguous, intent(in)  :: A(:, :)
                real(dp), contiguous, intent(out) :: B(:, :)
            endsubroutine uppertri_dp
        endinterface uppertri

        interface invert_permutation_old
            module subroutine invert_permutation_old(P, I)
                integer, intent(in)  :: P(:)
                integer, intent(out) :: I(size(P))
            endsubroutine invert_permutation_old
        endinterface invert_permutation_old

        interface permutecols
            module subroutine permutecols_dp(A, P, n)
                real(dp), intent(inout) :: A(n, n)
                integer , intent(inout) :: P(n)
                integer , intent(in)    :: n
            endsubroutine permutecols_dp
        endinterface permutecols

        interface permute_matrix_columns
            module subroutine permute_matrix_columns_dp(A, m, P, Pinv)
                real(dp), intent(inout) :: A(m, m)
                integer , intent(in)    :: m
                integer , intent(inout) :: P(m)
                integer , intent(inout) :: Pinv(m)
            endsubroutine permute_matrix_columns_dp
        endinterface permute_matrix_columns

        interface permute_matrix_rows
            module subroutine permute_matrix_rows_dp(A, m, P, Pinv)
                real(dp), intent(inout) :: A(m, m)
                integer , intent(in)    :: m
                integer , intent(inout) :: P(m)
                integer , intent(inout) :: Pinv(m)
            endsubroutine permute_matrix_rows_dp
        endinterface permute_matrix_rows

        interface rotate_matrix_columns
            module subroutine rotate_matrix_columns_dp(P, leader, A, m)
                integer , intent(inout) :: P(m)
                integer , intent(in)    :: leader
                real(dp), intent(inout) :: A(m, m)
                integer , intent(in)    :: m
            endsubroutine rotate_matrix_columns_dp
        endinterface rotate_matrix_columns

        interface rotate_matrix_rows
            module subroutine rotate_matrix_rows_dp(P, leader, A, m)
                integer , intent(inout) :: P(m)
                integer , intent(in)    :: leader
                real(dp), intent(inout) :: A(m, m)
                integer , intent(in)    :: m
            endsubroutine rotate_matrix_rows_dp
        endinterface rotate_matrix_rows

        interface permute
            module subroutine permute_dp(P, Pinv, x, m)
                integer , intent(inout) :: P(m)
                integer , intent(inout) :: Pinv(m)
                real(dp), intent(inout) :: x(m)
                integer , intent(in)    :: m
            endsubroutine permute_dp
        endinterface permute

        interface rotate
            module subroutine rotate_dp(P, leader, x)
                integer , intent(inout) :: P(:)
                integer , intent(in)    :: leader
                real(dp), intent(inout) :: x(:)
            endsubroutine rotate_dp
        endinterface rotate

        interface swap
            module subroutine swap_dp(x, i, j)
                real(dp), intent(inout) :: x(:)
                integer , intent(in)    :: i
                integer , intent(in)    :: j
            endsubroutine swap_dp
        endinterface swap

        interface invert
            module subroutine invert_dp(A, n, P, work, lwork, info, detsgn)
                real(dp), intent(inout)         :: A(n, n)
                integer , intent(in)            :: n
                integer , intent(out)           :: P(n)
                real(dp), intent(out)           :: work(lwork)
                integer , intent(in)            :: lwork
                integer , intent(out)           :: info
                integer , intent(out), optional :: detsgn
            endsubroutine invert_dp
        endinterface invert

        interface make_identity
            module subroutine make_identity_dp(A, n)
                real(dp), intent(inout) :: A(n, n)
                integer , intent(in)    :: n
            endsubroutine make_identity_dp
        endinterface make_identity

        interface zero_matrix
            module subroutine zero_matrix_dp(A, n)
                real(dp), intent(out) :: A(n, n)
                integer , intent(in)  :: n
            endsubroutine zero_matrix_dp
        endinterface zero_matrix

        interface copy_matrix
            module subroutine copy_matrix_dp(A, B)
                real(dp), intent(out) :: A(:, :)
                real(dp), intent(in)  :: B(:, :)
            endsubroutine copy_matrix_dp
        endinterface copy_matrix

        interface add_transpose
            module subroutine add_transpose_dp(A, B)
                real(dp), contiguous, intent(inout) :: A(:, :)
                real(dp), contiguous, intent(in)    :: B(:, :)
            endsubroutine add_transpose_dp
        endinterface add_transpose

        interface add_matrix
            module subroutine add_matrix_dp(A, B)
                real(dp), contiguous, intent(inout) :: A(:, :)
                real(dp), contiguous, intent(in)    :: B(:, :)
            endsubroutine add_matrix_dp
        endinterface add_matrix

        interface left_matmul
            module subroutine left_matmul_dp(A, B, work)
                real(dp), contiguous, intent(inout) :: A(:, :)
                real(dp), contiguous, intent(in)    :: B(:, :)
                real(dp), contiguous, intent(out)   :: work(:)
            endsubroutine left_matmul_dp
        endinterface left_matmul

        interface right_matmul
            module subroutine right_matmul_dp(A, B, work)
                real(dp), contiguous, intent(inout) :: A(:, :)
                real(dp), contiguous, intent(in)    :: B(:, :)
                real(dp), contiguous, intent(out)   :: work(:)
            endsubroutine right_matmul_dp
        endinterface right_matmul

        interface transpose
            module subroutine transpose_dp(A, B)
                real(dp), intent(out) :: A(:, :)
                real(dp), intent(in)  :: B(size(A, 2), size(A, 1))
            endsubroutine transpose_dp
        endinterface transpose

        interface laswpc
            module subroutine dlaswpc(n, a, lda, k1, k2, ipiv, incx)
                !    MODIFIED VERSION OF DLASWP
                !  -- LAPACK auxiliary routine --
                !  -- LAPACK is a software package provided by Univ. of Tennessee,    --
                !  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
                !
                !     .. Scalar Arguments ..
                integer            incx, k1, k2, lda, n
                !     ..
                !     .. Array Arguments ..
                integer            ipiv(*)
                double precision   A(lda, *)
            endsubroutine dlaswpc
        endinterface laswpc

        interface avgdiag
            module function avgdiag_dp(A, m) result(avg)
                real(dp), intent(in) :: A(m, m)
                integer , intent(in) :: m
                real(dp) :: avg
            endfunction avgdiag_dp
        endinterface avgdiag

        interface matdiff
            module function matdiff_dp(A, B, m, n, work, difftype) result(diff)
                real(dp)        , intent(in)            :: A(m, n)
                real(dp)        , intent(in)            :: B(m, n)
                integer         , intent(in)            :: m
                integer         , intent(in)            :: n
                real(dp)        , intent(out)           :: work(m, n)
                character(len=*), intent(in) , optional :: difftype
                real(dp) :: diff
            endfunction matdiff_dp
        endinterface matdiff

        interface determinant
            module subroutine determinant_dp(det, A, n, P, info)
                real(dp)                :: det
                real(dp), intent(inout) :: A(n, n)
                integer , intent(in)    :: n
                integer , intent(out)   :: P(n)
                integer , intent(out)   :: info
            endsubroutine determinant_dp
            module subroutine determinant_zp(det, A, n, P, info)
                complex(dp)                :: det
                complex(dp), intent(inout) :: A(n, n)
                integer    , intent(in)    :: n
                integer    , intent(out)   :: P(n)
                integer    , intent(out)   :: info
            endsubroutine determinant_zp
        endinterface determinant

        interface eigenvalues
            module subroutine eigenvalues_dp(A, wr, wi)
                real(dp), intent(in)  :: A(:, :)
                real(dp), intent(out) :: wr(size(A, 1))
                real(dp), intent(out) :: wi(size(A, 1))
            endsubroutine eigenvalues_dp
        endinterface eigenvalues

        interface diagonalize
            module subroutine diagonalize_dp(A, wr, wi, V)
                real(dp), intent(in)  :: A(:, :)
                real(dp), intent(out) :: wr(size(A, 1))
                real(dp), intent(out) :: wi(size(A, 1))
                real(dp), intent(out) :: V(size(A, 1), size(A, 1))
            endsubroutine diagonalize_dp
        endinterface diagonalize

        interface matrixswap
            module subroutine alloc_matrixswap_dp(A, B)
                real(dp), allocatable, intent(inout) :: A(:, :)
                real(dp), allocatable, intent(inout) :: B(:, :)
            endsubroutine alloc_matrixswap_dp
        endinterface matrixswap

        interface extractdiag
            module subroutine extractdiag_dp(A, D, m)
                real(dp), intent(in)  :: A(m, m)
                real(dp), intent(out) :: D(m)
                integer , intent(in)  :: m
            endsubroutine extractdiag_dp
        endinterface extractdiag

        interface id
            module subroutine id_dp(A, m)
                real(dp), intent(out) :: A(m, m)
                integer , intent(in)  :: m
            endsubroutine id_dp
        endinterface id

        interface qform
            module subroutine qform_dp(Q, tau, m, work, lwork)
                        real(dp), intent(inout) :: Q(m, m)
                        real(dp), intent(in)    :: tau(m)
                        integer , intent(in)    :: m
                        real(dp), intent(out)   :: work(lwork)
                        integer , intent(in)    :: lwork
            endsubroutine qform_dp
        endinterface qform

        interface append_column
            module subroutine append_column_dp(A, x)
                real(dp), allocatable, intent(inout) :: A(:, :)
                real(dp),              intent(in)    :: x(:)
            endsubroutine append_column_dp
        endinterface append_column

        interface lu
            module subroutine lu_dp(A, piv)
                real(dp), contiguous, intent(inout) :: A(:, :)
                integer,  contiguous, intent(out)   :: piv(:)
            endsubroutine lu_dp
        endinterface lu

        interface pivsgn
            module function pivsgn(piv, pivl) result(sgn)
                integer, intent(in)           :: piv(:)
                integer, intent(in), optional :: pivl
                integer :: sgn
            endfunction pivsgn
        endinterface pivsgn

        interface ludet
            module function ludet_dp(A, piv) result(det)
                real(dp), intent(in) :: A(:, :)
                integer,  intent(in) :: piv(:)
                real(dp) :: det
            endfunction ludet_dp
        endinterface ludet

        interface luinv
            module subroutine luinv_dp(A, piv, work)
                real(dp), contiguous, intent(inout) :: A(:, :)
                integer , contiguous, intent(in)    :: piv(:)
                real(dp), contiguous, intent(out)   :: work(:)
            endsubroutine luinv_dp
        endinterface luinv




    ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    ! ||||||||||BLAS and LAPACK interfaces ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    ! -----------------------------------------------------------------------------------------------------
    ! Scalar operations -----------------------------------------------------------------------------------

    interface
        real(sp) function scabs1(z)
            import                  :: sp
            complex(sp), intent(in) :: z
        endfunction scabs1
    endinterface

    interface
        real(dp) function dcabs1(z)
            import                  :: dp
            complex(dp), intent(in) :: z
        endfunction dcabs1
    endinterface

    ! end Scalar operations -------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Level 1 BLAS: vector ops ----------------------------------------------------------------------------
    
    interface
        subroutine dcopy(n, dx, incx, dy, incy)
            import                 :: dp
            integer , intent(in)   :: n
            real(dp), intent(in)   :: dx(*)
            integer , intent(in)   :: incx
            real(dp), intent(out)  :: dy(*)
            integer , intent(in)   :: incy
        endsubroutine dcopy
    endinterface

    interface
        subroutine dscal(n, da, dx, incx)
            import                  :: dp
            integer , intent(in)    :: n
            real(dp), intent(in)    :: da
            real(dp), intent(inout) :: dx(*)
            integer , intent(in)    :: incx
        endsubroutine dscal
    endinterface

    interface
        subroutine dswap(n, dx, incx, dy, incy)
            import                  :: dp
            integer , intent(in)    :: n
            real(dp), intent(inout) :: dx(*)
            integer , intent(in)    :: incx
            real(dp), intent(inout) :: dy(*)
            integer , intent(in)    :: incy
        endsubroutine dswap
    endinterface

    interface
        subroutine daxpy(n   , da, dx, incx, dy,   &
                         incy)
            import                  :: dp
            integer , intent(in)    :: n
            real(dp), intent(in)    :: da
            real(dp), intent(in)    :: dx(*)
            integer , intent(in)    :: incx
            real(dp), intent(inout) :: dy(*)
            integer , intent(in)    :: incy
        endsubroutine daxpy
    endinterface

    ! end Level 1 BLAS: vector ops ------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Level 2 BLAS: matrix-vector ops ---------------------------------------------------------------------
    
    interface
        subroutine dger(m, n   , alpha, x  , incx,  &
                        y, incy, a    , lda)
            import                  :: dp
            integer , intent(in)    :: m
            integer , intent(in)    :: n
            real(dp), intent(in)    :: alpha
            real(dp), intent(in)    :: x(*)
            integer , intent(in)    :: incx
            real(dp), intent(in)    :: y(*)
            integer , intent(in)    :: incy
            real(dp), intent(inout) :: a(lda, *)
            integer , intent(in)    :: lda
        endsubroutine dger
    endinterface

    ! end Level 2 BLAS: matrix-vector ops -----------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Level 3 BLAS: matrix-matrix ops ---------------------------------------------------------------------
    
    interface
        subroutine dgemm(transa, transb, m  , n, k  ,   &
                         alpha , a     , lda, b, ldb,   &
                         beta  , c     , ldc)
            import                   :: dp
            character, intent(in)    :: transa
            character, intent(in)    :: transb
            integer  , intent(in)    :: m
            integer  , intent(in)    :: n
            integer  , intent(in)    :: k
            real(dp) , intent(in)    :: alpha
            real(dp) , intent(in)    :: a(lda, *)
            integer  , intent(in)    :: lda
            real(dp) , intent(in)    :: b(ldb, *)
            integer  , intent(in)    :: ldb
            real(dp) , intent(in)    :: beta
            real(dp) , intent(inout) :: c(ldc, *)
            integer  , intent(in)    :: ldc
        endsubroutine dgemm
    endinterface

    interface
        subroutine dtrmm(side, uplo , transa, diag, m,      &
                         n   , alpha, a     , lda , b,      &
                         ldb)
            import                   :: dp
            character, intent(in)    :: side
            character, intent(in)    :: uplo
            character, intent(in)    :: transa
            character, intent(in)    :: diag
            integer  , intent(in)    :: m
            integer  , intent(in)    :: n
            real(dp) , intent(in)    :: alpha
            real(dp) , intent(in)    :: a(lda, *)
            integer  , intent(in)    :: lda
            real(dp) , intent(inout) :: b(ldb, *)
            integer  , intent(in)    :: ldb
        endsubroutine dtrmm
    endinterface

    ! end Level 3 BLAS: matrix-matrix ops -----------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Linear solve, AX = B --------------------------------------------------------------------------------

    interface
        subroutine dgetrf(m, n, a, lda, ipiv, info)
            import                  :: dp
            integer , intent(in)    :: m
            integer , intent(in)    :: n
            real(dp), intent(inout) :: a(lda, *)
            integer , intent(in)    :: lda
            integer , intent(out)   :: ipiv(*)
            integer , intent(out)   :: info
        endsubroutine dgetrf
    endinterface

    interface
        subroutine zgetrf(m   , n, a, lda, ipiv,    &
                          info)
            import                     :: dp
            integer    , intent(in)    :: m
            integer    , intent(in)    :: n
            complex(dp), intent(inout) :: a(lda, *)
            integer    , intent(in)    :: lda
            integer    , intent(out)   :: ipiv(*)
            integer    , intent(out)   :: info
        endsubroutine zgetrf
    endinterface

    interface
        subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
            import                  :: dp
            integer , intent(in)    :: n
            real(dp), intent(inout) :: a(lda, *)
            integer , intent(in)    :: lda
            integer , intent(in)    :: ipiv(*)
            real(dp), intent(out)   :: work(*)
            integer , intent(in)    :: lwork
            integer , intent(out)   :: info
        endsubroutine dgetri
    endinterface

    ! end Linear solve, AX = B ----------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Least squares ---------------------------------------------------------------------------------------

    ! end Least squares -----------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Orthogonal/unitary factors (QR, CS, etc.) -----------------------------------------------------------

    interface
        subroutine dgeqp3(m  , n   , a    , lda , jpvt,     &
                          tau, work, lwork, info)
            import                  :: dp
            integer , intent(in)    :: m
            integer , intent(in)    :: n
            real(dp), intent(inout) :: A(lda, *)
            integer , intent(in)    :: lda
            integer , intent(inout) :: jpvt(*)
            real(dp), intent(out)   :: tau(*)
            real(dp), intent(out)   :: work(*)
            integer , intent(in)    :: lwork
            integer , intent(out)   :: info
        endsubroutine dgeqp3
    endinterface

    interface
        subroutine dorgqr(m  , n   , k    , a   , lda,      &
                          tau, work, lwork, info)
            import                  :: dp
            integer , intent(in)    :: m
            integer , intent(in)    :: n
            integer , intent(in)    :: k
            real(dp), intent(inout) :: A(lda, *)
            integer , intent(in)    :: lda
            real(dp), intent(in)    :: tau(*)
            real(dp), intent(out)   :: work(*)
            integer , intent(in)    :: lwork
            integer , intent(out)   :: info
        endsubroutine dorgqr
    endinterface

    interface
        subroutine dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
            import                   :: dp
            character, intent(in)    :: side
            character, intent(in)    :: trans
            integer  , intent(in)    :: m
            integer  , intent(in)    :: n
            integer  , intent(in)    :: k
            real(dp) , intent(in)    :: a(lda, *)
            integer  , intent(in)    :: lda
            real(dp) , intent(in)    :: tau(*)
            real(dp) , intent(inout) :: c(ldc, *)
            integer  , intent(in)    :: ldc
            real(dp) , intent(in)    :: work(*)
            integer  , intent(in)    :: lwork
            integer  , intent(out)   :: info
        endsubroutine dormqr
    endinterface

    ! end Orthogonal/unitary factors (QR, CS, etc.) -------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Non-symmetric eigenvalues ---------------------------------------------------------------------------

    interface
        subroutine dgeev(jobvl, jobvr, n    , a   , lda,    &
                         wr   , wi   , vl   , ldvl, vr ,    &
                         ldvr , work , lwork, info)
            import                   :: dp
            character, intent(in)    :: jobvl
            character, intent(in)    :: jobvr
            integer  , intent(in)    :: n
            real(dp) , intent(inout) :: a(lda, *)
            integer  , intent(in)    :: lda
            real(dp) , intent(out)   :: wr(*)
            real(dp) , intent(out)   :: wi(*)
            real(dp) , intent(out)   :: vl(ldvl, *)
            integer  , intent(in)    :: ldvl
            real(dp) , intent(out)   :: vr(ldvr, *)
            integer  , intent(in)    :: ldvr
            real(dp) , intent(out)   :: work(*)
            integer  , intent(in)    :: lwork
            integer  , intent(out)   :: info
        endsubroutine dgeev
    endinterface

    ! end Non-symmetric eigenvalues -----------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Hermitian/symmetric eigenvalues ---------------------------------------------------------------------

    ! end Hermitian/symmetric eigenvalues -----------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Singular Value Decomposition (SVD) ------------------------------------------------------------------
    interface
        subroutine dgesdd(jobz, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, iwork, info)
            import                   :: dp
            character, intent(in)    :: jobz
            integer  , intent(in)    :: m
            integer  , intent(in)    :: n
            real(dp) , intent(inout) :: a(lda, *)
            integer  , intent(in)    :: lda
            real(dp) , intent(out)   :: s(*)
            real(dp) , intent(out)   :: u(ldu, *)
            integer  , intent(in)    :: ldu
            real(dp) , intent(out)   :: vt(ldvt, *)
            integer  , intent(in)    :: ldvt
            real(dp) , intent(out)   :: work(*)
            integer  , intent(in)    :: lwork
            integer  , intent(out)   :: iwork(*)
            integer  , intent(out)   :: info
        endsubroutine dgesdd
    endinterface

    interface
        subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
            import                   :: dp
            character, intent(in)    :: jobu
            character, intent(in)    :: jobvt
            integer  , intent(in)    :: m
            integer  , intent(in)    :: n
            real(dp) , intent(inout) :: a(lda, *)
            integer  , intent(in)    :: lda
            real(dp) , intent(out)   :: s(*)
            real(dp) , intent(out)   :: u(ldu, *)
            integer  , intent(in)    :: ldu
            real(dp) , intent(out)   :: vt(ldvt, *)
            integer  , intent(in)    :: ldvt
            real(dp) , intent(out)   :: work(*)
            integer  , intent(in)    :: lwork
            integer  , intent(out)   :: info
        endsubroutine dgesvd
    endinterface

    interface
        subroutine dgesvj(joba, jobu, jobv, m, n, a, lda, sva, mv, v, ldv, work, lwork, info)
            import :: dp
            character, intent(in)    :: joba
            character, intent(in)    :: jobu
            character, intent(in)    :: jobv
            integer  , intent(in)    :: m
            integer  , intent(in)    :: n
            real(dp) , intent(inout) :: a(lda, *)
            integer  , intent(in)    :: lda
            real(dp) , intent(out)   :: sva(*)
            integer  , intent(in)    :: mv
            real(dp) , intent(inout) :: v(ldv, *)
            integer  , intent(in)    :: ldv
            real(dp) , intent(inout) :: work(*)
            integer  , intent(in)    :: lwork
            integer  , intent(out)   :: info
        endsubroutine dgesvj
    endinterface
    ! end Singular Value Decomposition (SVD) --------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! BLAS-like -------------------------------------------------------------------------------------------

    interface
        subroutine dlaset(uplo, m  , n, alpha, beta,  &
                          a   , lda)
            import                 :: dp
            character, intent(in)  :: uplo
            integer  , intent(in)  :: m
            integer  , intent(in)  :: n
            real(dp) , intent(in)  :: alpha
            real(dp) , intent(in)  :: beta
            real(dp) , intent(out) :: a(lda, *)
            integer  , intent(in)  :: lda
        endsubroutine dlaset
    endinterface

    interface
        subroutine dlascl2(m, n, d, x, ldx)
            import                  :: dp
            integer , intent(in)    :: m
            integer , intent(in)    :: n
            real(dp), intent(in)    :: d(*)
            real(dp), intent(inout) :: x(ldx, *)
            integer , intent(in)    :: ldx
        endsubroutine dlascl2
    endinterface

    interface
        subroutine dlacpy(uplo, m  , n, a, lda,   &
                          b   , ldb)
            import                 :: dp
            character, intent(in)  :: uplo
            integer  , intent(in)  :: m
            integer  , intent(in)  :: n
            real(dp) , intent(in)  :: a(lda, *)
            integer  , intent(in)  :: lda
            real(dp) , intent(out) :: b(ldb, *)
            integer  , intent(in)  :: ldb
        endsubroutine dlacpy
    endinterface

    interface
        real(dp) function dlange(norm, m, n, a, lda,    &
                                 work)
            import                 :: dp
            character, intent(in)  :: norm
            integer  , intent(in)  :: m
            integer  , intent(in)  :: n
            real(dp) , intent(in)  :: a(lda, *)
            integer  , intent(in)  :: lda
            real(dp) , intent(out) :: work(*)
        endfunction dlange
    endinterface

    ! end BLAS-like ---------------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------


    ! -----------------------------------------------------------------------------------------------------
    ! Auxiliary routines ----------------------------------------------------------------------------------

    ! end Auxiliary routines ------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------
endmodule stdlinalg