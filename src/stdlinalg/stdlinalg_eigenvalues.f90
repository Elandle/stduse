submodule(stdlinalg) stdlinalg_eigenvalues
    implicit none
    contains
        module procedure eigenvalues_dp
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
        endprocedure eigenvalues_dp

        module procedure diagonalize_dp
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
        endprocedure diagonalize_dp






endsubmodule stdlinalg_eigenvalues