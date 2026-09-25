module stdlinalg
    use stduse
    implicit none

    interface alloc_matrixswap
        module procedure alloc_matrixswap_dp
    endinterface alloc_matrixswap

    interface extractdiag
        module procedure extractdiag_dp
    endinterface

    interface qform
        module procedure qform_dp
    endinterface

    interface id
        module procedure id_dp
    endinterface id


    interface append_column
        module procedure :: append_column_dp
    endinterface append_column


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




























    ! -----------------------------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------
    ! BLAS and LAPACK interfaces --------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------------------------------

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
        subroutine dgetrf(m   , n, a, lda, ipiv,   &
                          info)
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
        subroutine dgetri(n    , a   , lda, ipiv, work,    &
                          lwork, info)
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
            integer , intent(in)    :: d(*)
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































    contains

        real(dp) function twonorm(A) result(norm)
            !
            ! Returns the 2 norm of the matrix A
            !
            real(dp), intent(in)  :: A(:, :)
            real(dp), allocatable :: B(:, :)
            real(dp), allocatable :: S(:)
            real(dp), allocatable :: work(:)
            integer               :: m
            integer               :: n
            integer               :: lwork
            integer               :: info
            
            m = size(A, 1)
            n = size(A, 2)
            lwork = 5 * (3 * min(m, n) + max(m, n) + 5 * min(m, n))

            allocate(B(m, n))
            allocate(S(min(m, n)))
            allocate(work(lwork))

            call dlacpy('A', m, n, A, m, B, m)
            call dgesvd('N', 'N', m, n, B, m, S, work, lwork, work, lwork, work, lwork, info)

            norm = S(1) 
        endfunction twonorm



        !> \brief Updates \f$A = AD\f$ for a square matrix \f$A\f$ and diagonal matrix \f$D\f$ stored as a vector.
        !!
        !! \param[inout] A (`real(dp), dimension(n, n)`) \f$n\times n\f$ matrix \f$A\f$ to update.
        !! \param[in]    D (`real(dp), dimension(n)`)    Diagonal matrix \f$D\f$ stored as a vector.
        !! \param[in]    n (`integer`)                   Dimension of \f$A\f$ and \f$D\f$.
        !! \see Todo: maybe change the `do` to a `do concurrent`?
        subroutine right_diagmult(A, D)
            !
            ! Updates:
            !
            ! A = A * diag(D)
            !
            ! where A is an m x n matrix and D is an n-long vector.
            !
            real(dp), contiguous, intent(inout) :: A(:, :)
            real(dp), contiguous, intent(in)    :: D(:)

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
        endsubroutine right_diagmult

        !> \brief Updates \f$A = DA\f$ for a square matrix \f$A\f$ and diagonal matrix \f$D\f$ stored as a vector.
        !!
        !! \param[inout] A (`real(dp), dimension(n, n)`) \f$n\times n\f$ matrix \f$A\f$ to update.
        !! \param[in]    D (`real(dp), dimension(n)`)    Diagonal matrix \f$D\f$ stored as a vector.
        !! \param[in]    n (`integer`)                   Dimension of \f$A\f$ and \f$D\f$.
        !! \see Todo: maybe change the `do` to a `do concurrent`? </p>
        !! Some BLAS/LAPACK distributions contain the `dlascl2` subroutine and some do not.
        !! `dlascl2` should be used if possible. If not, provided alternative code can be used.
        subroutine left_diagmult(A, D)
            !
            ! Updates:
            !
            ! A = diag(D) * A
            !
            ! where A is an m x n matrix and D is an m-long vector.
            !
            real(dp), contiguous, intent(inout) :: A(:, :)
            real(dp), contiguous, intent(in)    :: D(:)

            call left_diagmult_helper(A, D, size(A, 1), size(A, 2))
        contains
            subroutine left_diagmult_helper(A, D, m, n)
                integer , intent(in)    :: m
                integer , intent(in)    :: n
                real(dp), intent(inout) :: A(m, n)
                real(dp), intent(in)    :: D(m)

                integer :: i

                ! 
                ! call dlascl2(m, n, D, A, m)

                ! Scale row i of A by D(i)
                do i = 1, m
                    call dscal(n, D(i), A(i, 1), m)
                enddo
            endsubroutine left_diagmult_helper
        endsubroutine left_diagmult

        !> \brief Updates \f$A = AD^{-1}\f$ for a square matrix \f$A\f$ and diagonal matrix \f$D\f$ stored as a vector.
        !!
        !! \param[inout] A (`real(dp), dimension(n, n)`) \f$n\times n\f$ matrix \f$A\f$ to update.
        !! \param[in]    D (`real(dp), dimension(n)`)    Diagonal matrix \f$D\f$ stored as a vector.
        !! \param[in]    n (`integer`)                   Dimension of \f$A\f$ and \f$D\f$.
        !! \see Todo: maybe change the `do` to a `do concurrent`? </p>
        !! This code does not check if \f$D\f$ has all nonzero diagonal entries (in which case
        !! \f$D^{-1}\f$ does not exist).
        subroutine right_diaginvmult(A, D)
            !
            ! Updates:
            !
            ! A = A * inv(diag(D))
            !
            ! where A is an m x n matrix and D is an n-long vector.
            !
            real(dp), contiguous, intent(inout) :: A(:, :)
            real(dp), contiguous, intent(in)    :: D(:)

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
        endsubroutine right_diaginvmult

        !> \brief Updates \f$A = D^{-1}A\f$ for a square matrix \f$A\f$ and diagonal matrix \f$D\f$ stored as a vector.
        !!
        !! \param[inout] A (`real(dp), dimension(n, n)`) \f$n\times n\f$ matrix \f$A\f$ to update.
        !! \param[in]    D (`real(dp), dimension(n)`)    Diagonal matrix \f$D\f$ stored as a vector.
        !! \param[in]    n (`integer`)                   Dimension of \f$A\f$ and \f$D\f$.
        !! \see Todo: maybe change the `do` to a `do concurrent`? </p>
        !! This code does not check if \f$D\f$ has all nonzero diagonal entries (in which case
        !! \f$D^{-1}\f$ does not exist).
        !! Some BLAS/LAPACK distributions contain the `dlascl2` subroutine and some do not.
        !! `dlascl2` should be used if possible. If not, provided alternative code can be used.
        subroutine left_diaginvmult(A, D)
            !
            ! Updates:
            !
            ! A = inv(diag(D)) * A
            !
            ! where A is an m x n matrix and D is an m-vector.
            !
            real(dp), contiguous, intent(inout) :: A(:, :)
            real(dp), contiguous, intent(in)    :: D(:)

            call left_diaginvmult_helper(A, D, size(A, 1), size(A, 2))
        contains
            subroutine left_diaginvmult_helper(A, D, m, n)
                integer , intent(in)    :: m
                integer , intent(in)    :: n
                real(dp), intent(inout) :: A(m, n)
                real(dp), intent(in)    :: D(m)

                integer :: i

                ! call dlarscl2(m, n, D, A, m)

                ! Scale row i of A by 1/D(i)
                do i = 1, m
                    call dscal(n, 1.0_dp / D(i), A(i, 1), m)
                enddo
            endsubroutine left_diaginvmult_helper
        endsubroutine left_diaginvmult

        !> \brief Sets \f$D = \text{diag}(A)\f$ for a square matrix \f$A\f$ and diagonal matrix \f$D\f$ stored as a vector.
        !!
        !! \param[inout] A (`real(dp), dimension(n, n)`) \f$n\times n\f$ matrix \f$A\f$ to update.
        !! \param[in]    D (`real(dp), dimension(n)`)    Diagonal matrix \f$D\f$ stored as a vector.
        !! \param[in]    n (`integer`)                   Dimension of \f$A\f$ and \f$D\f$.
        subroutine diag(A, D)
            !
            ! Sets:
            !
            ! D = diagonal(A)
            !
            ! where A is an m x n matrix and D has length min(m, n).
            !
            real(dp), contiguous, intent(in)  :: A(:, :)
            real(dp), contiguous, intent(out) :: D(:)

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
        endsubroutine diag

        subroutine uppertri(A, B)
            !
            ! Sets:
            !
            ! B = uppertri(A)
            !
            ! That is, B is set to be the upper triangular part of A (including the diagonal),
            ! with all other entries set to zero.
            !
            real(dp), contiguous, intent(in)  :: A(:, :)
            real(dp), contiguous, intent(out) :: B(:, :)

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
        endsubroutine uppertri

        
        subroutine invert_permutation_old(P, I, n)
            !
            ! Sets:
            !
            ! I = inv(P)
            !
            ! In terms of permutation matrices stored as length n vectors
            !
            integer, intent(in)  :: P(n)
            integer, intent(out) :: I(n)
            integer, intent(in)  :: n

            integer :: j

            do j = 1, n
                I(P(j)) = j ! Better way to do this using array indexing?
            enddo

            ! Maybe:
            ! I(P) = [(j, j = 1, n)]
        endsubroutine invert_permutation_old
        
        subroutine permutecols(A, P, n)
            !
            ! Updates:
            !
            ! A = A * P
            !
            ! Where P is a permutation matrix stored as a vector.
            !
            ! This subroutine is implemented basically as a copy of Julia's intrinsic permutecols
            !
            real(dp), intent(inout) :: A(n, n)
            integer , intent(inout) :: P(n)
            integer , intent(in)    :: n

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
        endsubroutine permutecols

        ! These following permutation algorithms are all based off the same underlying algorithm
        ! TODO:
        ! Find an inplace permutation algorithm that doesn't need an inverse permutation

        subroutine permute_matrix_columns(A, m, P, Pinv)
            real(dp), intent(inout) :: A(m, m)
            integer , intent(in)    :: m
            integer , intent(inout) :: P(m)
            integer , intent(inout) :: Pinv(m)

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
        endsubroutine permute_matrix_columns

        subroutine permute_matrix_rows(A, m, P, Pinv)
            real(dp), intent(inout) :: A(m, m)
            integer , intent(in)    :: m
            integer , intent(inout) :: P(m)
            integer , intent(inout) :: Pinv(m)

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
        endsubroutine permute_matrix_rows

        subroutine rotate_matrix_columns(P, leader, A, m)
            ! Part of permute
            integer , intent(inout) :: P(m)
            integer , intent(in)    :: leader
            real(dp), intent(inout) :: A(m, m)
            integer , intent(in)    :: m

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
        endsubroutine rotate_matrix_columns

        subroutine rotate_matrix_rows(P, leader, A, m)
            ! Part of permute
            integer , intent(inout) :: P(m)
            integer , intent(in)    :: leader
            real(dp), intent(inout) :: A(m, m)
            integer , intent(in)    :: m

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
        endsubroutine rotate_matrix_rows

        subroutine permute(P, Pinv, x, m)
            ! Fich, Munro, and Poblete
            ! Permuting in Place
            ! Figure 5
            integer , intent(inout) :: P(m)
            integer , intent(inout) :: Pinv(m)
            real(dp), intent(inout) :: x(m)
            integer , intent(in)    :: m

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
        endsubroutine permute

        subroutine rotate(P, leader, x)
            ! Part of permute
            integer , intent(inout) :: P(:)
            integer , intent(in)    :: leader
            real(dp), intent(inout) :: x(:)

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
        endsubroutine rotate

        subroutine swap(x, i, j)
            real(dp), intent(inout) :: x(:)
            integer , intent(in)    :: i
            integer , intent(in)    :: j

            real(dp) :: temp

            temp = x(i)
            x(i) = x(j)
            x(j) = temp
        endsubroutine swap

        subroutine invert(A, n, P, work, lwork, info, detsgn)
            !
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
            real(dp), intent(inout)         :: A(n, n)
            integer , intent(in)            :: n
            integer , intent(out)           :: P(n)
            real(dp), intent(out)           :: work(lwork)
            integer , intent(in)            :: lwork
            integer , intent(out)           :: info
            integer , intent(out), optional :: detsgn

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
        endsubroutine invert

        subroutine make_identity(A, n)
            !
            ! Sets:
            !
            ! A = id
            !
            ! where id is the identity matrix and A is an n x n matrix
            !
            real(dp), intent(inout) :: A(n, n)
            integer , intent(in)    :: n

            ! A = id
            call dlaset('A', n, n, 0.0_dp, 1.0_dp, A, n)
        endsubroutine make_identity

        subroutine zero_matrix(A, n)
            !
            ! Sets:
            !
            ! A = 0
            !
            real(dp), intent(out) :: A(n, n)
            integer , intent(in)  :: n

            ! A = 0
            call dlaset('A', n, n, 0.0_dp, 0.0_dp, A, n)
        endsubroutine zero_matrix

        !> Copies \f$A = B\f$, where \f$A\f$ and \f$B\f$ are both \f$m\times n\f$ matrices.
        !!
        !! \param[out]  A  (`real(dp), dimension(:, :)`) Matrix to copy into.
        !! \param[in]   B  (`real(dp), dimension(size(A, 1), size(A, 2))`) Matrix to copy from.
        subroutine copy_matrix(A, B)
            real(dp), intent(out) :: A(:, :)
            real(dp), intent(in)  :: B(:, :)

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
        endsubroutine copy_matrix

        subroutine add_transpose(A, B)
            !
            ! Updates:
            !
            ! A = A + transpose(B)
            !
            ! where A is m x n and B is n x m.
            !
            real(dp), contiguous, intent(inout) :: A(:, :)
            real(dp), contiguous, intent(in)    :: B(:, :)

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
        end subroutine add_transpose

        subroutine add_matrix(A, B)
            !
            ! Updates:
            !
            ! A = A + B
            !
            ! where A and B are m x n matrices
            !
            real(dp), contiguous, intent(inout) :: A(:, :)
            real(dp), contiguous, intent(in)    :: B(:, :)

            call add_matrix_helper(A, B, size(A, 1), size(A, 2))

            contains
                subroutine add_matrix_helper(A, B, m, n)
                    integer , intent(in)    :: m
                    integer , intent(in)    :: n
                    real(dp), intent(inout) :: A(m, n)
                    real(dp), intent(in)    :: B(m, n)

                    call daxpy(size(A, 1) * size(A, 2), 1.0_dp, B, 1, A, 1)
                endsubroutine add_matrix_helper
        endsubroutine add_matrix

        subroutine left_matmul(A, B, work)
            !
            ! Updates:
            !
            ! A = B * A
            !
            ! where A is an m x n matrix and B an m x m matrix (so the result can still be stored in A).
            ! work is a workspace array (any dimensional, eg 1 or 2 as long as it has enough space)
            ! at least m * n long.
            !
            real(dp), contiguous, intent(inout) :: A(:, :)
            real(dp), contiguous, intent(in)    :: B(:, :)
            real(dp), contiguous, intent(out)   :: work(:)

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
        endsubroutine left_matmul

        subroutine right_matmul(A, B, work)
            !
            ! Updates:
            !
            ! A = A * B
            !
            ! where A is an m x n matrix and B an n x n matrix (so the result can still be stored in A).
            ! work is a workspace array (any dimensional, eg 1 or 2 as long as it has enough space)
            ! at least m * n long.
            !
            real(dp), contiguous, intent(inout) :: A(:, :)
            real(dp), contiguous, intent(in)    :: B(:, :)
            real(dp), contiguous, intent(out)   :: work(:)

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
        endsubroutine right_matmul

        subroutine transpose(A, B)
            !
            ! Sets A to the transpose of B.
            !
            real(dp), intent(out) :: A(:, :)
            real(dp), intent(in)  :: B(size(A, 2), size(A, 1))

            integer :: i, j, m, n

            m = size(A, 1) ; n = size(A, 2)

            ! I believe Intel's MKL has something for this.

            do j = 1, n
                do i = 1, n
                    A(i, j) = B(j, i)
                enddo
            enddo
        endsubroutine transpose

        subroutine dlaswpc(n, a, lda, k1, k2, ipiv, incx)
            !
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
            !     ..
            !
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
        endsubroutine dlaswpc

        real(dp) function avgdiag(A, m) result(avg)
            real(dp), intent(in) :: A(m, m)
            integer , intent(in) :: m

            integer  :: j

            avg = 0.0_dp
            do j = 1, m
                avg = avg + A(j, j)
            enddo
            avg = avg / m
        endfunction avgdiag

        function matdiff(A, B, m, n, work, difftype) result(diff)
            real(dp)        , intent(in)            :: A(m, n)
            real(dp)        , intent(in)            :: B(m, n)
            integer         , intent(in)            :: m
            integer         , intent(in)            :: n
            real(dp)        , intent(out)           :: work(m, n)
            character(len=*), intent(in) , optional :: difftype

            ! real(dp)    , external    :: dlange
            real(dp)                  :: diff
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
        endfunction matdiff

        subroutine determinant(det, A, n, P, info)
            !
            ! Sets:
            !
            ! det = det(A)
            !
            ! by first factorising A = P * L * U (destroying A),
            ! and calculating the determinant by multiplying the diagonal of U
            ! and adjusting the sign based on P
            !
            real(dp)                :: det
            real(dp), intent(inout) :: A(n, n)
            integer , intent(in)    :: n
            integer , intent(out)   :: P(n)
            integer , intent(out)   :: info

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
        endsubroutine determinant

        subroutine zdeterminant(det, A, n, P, info)
            !
            ! Sets:
            !
            ! det = det(A)
            !
            ! by first factorising A = P * L * U (destroying A),
            ! and calculating the determinant by multiplying the diagonal of U
            ! and adjusting the sign based on P
            !
            complex(dp)                :: det
            complex(dp), intent(inout) :: A(n, n)
            integer    , intent(in)    :: n
            integer    , intent(out)   :: P(n)
            integer    , intent(out)   :: info

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
        endsubroutine zdeterminant

        subroutine eigenvalues(A, wr, wi)
            real(dp), intent(in)  :: A(:, :)
            real(dp), intent(out) :: wr(size(A, 1))
            real(dp), intent(out) :: wi(size(A, 1))

            integer :: m
            real(dp), allocatable :: B(:, :)
            real(dp), allocatable :: work(:)
            integer               :: lwork
            integer               :: info

            m = size(A, 1)
            lwork = 8 * m
            allocate(B(m, m))
            allocate(work(lwork))
            call copy_matrix(B, A)


            call dgeev('N', 'N', m, B, m, wr, wi, B, m, B, m, work, lwork, info)
        endsubroutine eigenvalues


        subroutine diagonalize(A, wr, wi, V)
            real(dp), intent(in)  :: A(:, :)
            real(dp), intent(out) :: wr(size(A, 1))
            real(dp), intent(out) :: wi(size(A, 1))
            real(dp), intent(out) :: V(size(A, 1), size(A, 1))

            integer :: m
            real(dp), allocatable :: B(:, :)
            real(dp), allocatable :: work(:)
            integer               :: lwork
            integer               :: info

            m = size(A, 1)
            lwork = 8 * m
            allocate(B(m, m))
            allocate(work(lwork))
            call copy_matrix(B, A)


            call dgeev('N', 'V', m, B, m, wr, wi, B, m, V, m, work, lwork, info)
        endsubroutine diagonalize


        subroutine alloc_matrixswap_dp(A, B)
            !
            ! Effectively interchanges the names of the allocatable real(dp) matrices A and B.
            ! If A is matrix m1 and B matrix m2, then after calling this
            ! A is matrix m2 and B is matrix m1.
            !
            real(dp), allocatable, intent(inout) :: A(:, :)
            real(dp), allocatable, intent(inout) :: B(:, :)

            real(dp), allocatable :: temp(:, :)

            call move_alloc(A, temp)
            call move_alloc(B, A)
            call move_alloc(temp, B)
        endsubroutine alloc_matrixswap_dp

        subroutine extractdiag_dp(A, D, m)
            real(dp), intent(in)  :: A(m, m)
            real(dp), intent(out) :: D(m)
            integer , intent(in)  :: m

            integer :: i

            do i = 1, m
                D(i) = A(i, i)
            enddo
        endsubroutine extractdiag_dp

        subroutine id_dp(A, m)
            real(dp), intent(out) :: A(m, m)
            integer , intent(in)  :: m

            call dlaset('A', m, m, 0.0_dp, 1.0_dp, A, m)
        endsubroutine id_dp

        subroutine qform_dp(Q, tau, m, work, lwork)
            !
            ! Wrapper for calling LAPACK's real(dp) Q formation routine dorgqr.
            ! Forms the matrix Q in a QR or QRP factorization.
            ! The information about Q should be contained in Q, tau on input,
            ! and on output Q is changed to contain the full Q.
            ! 
            real(dp), intent(inout) :: Q(m, m)
            real(dp), intent(in)    :: tau(m)
            integer , intent(in)    :: m
            real(dp), intent(out)   :: work(lwork)
            integer , intent(in)    :: lwork

            integer :: info

            call dorgqr(m, m, m, Q, m, tau, work, lwork, info)
        endsubroutine qform_dp

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











endmodule stdlinalg