submodule(stduse) stduse_print
    implicit none
    ! Default format statements for printing.
    character(len=*), parameter :: dmatrixfmt = "(f17.8)"
    character(len=*), parameter :: dvectorfmt = "(f17.8)"
    character(len=*), parameter :: ivectorfmt = "(i6)"
    character(len=*), parameter :: imatrixfmt = "(i6)"
    contains
        module procedure print_dmatrix
            integer                       :: m, n, i, j
            character(len=:), allocatable :: afmt

            if (present(fmt)) then
                afmt = fmt
            else
                afmt = dmatrixfmt
            endif

            m = size(A, 1) ; n = size(A, 2)

            if (present(message)) write (ounit, "(a)") message

            do i = 1, m
                do j = 1, n
                    write(unit=ounit, fmt=afmt, advance="no") A(i, j)
                enddo
                write(unit=ounit, fmt="(a)") ""
            enddo
        endprocedure print_dmatrix

        module procedure print_imatrix
            integer                       :: m, n, i, j
            character(len=:), allocatable :: afmt

            if (present(fmt)) then
                afmt = fmt
            else
                afmt = imatrixfmt
            endif

            m = size(A, 1) ; n = size(A, 2)

            if (present(message)) write (ounit, "(a)") message

            do i = 1, m
                do j = 1, n
                    write(unit=ounit, fmt=afmt, advance="no") A(i, j)
                enddo
                write(unit=ounit, fmt="(a)") ""
            enddo
        endprocedure print_imatrix

        module procedure print_dvector
            integer                       :: m, i
            character(len=:), allocatable :: afmt
            logical                       :: aadvance

            if (present(advance)) then
                aadvance = advance
            else
                aadvance = .false.
            endif

            if (present(fmt)) then
                afmt = fmt
            else
                afmt = dvectorfmt
            endif

            m = size(v)

            if (present(message)) write (ounit, "(a)") message

            if (aadvance) then
                do i = 1, m
                    write(unit=ounit, fmt=afmt) v(i)
                enddo
            else
                do i = 1, m
                    write(unit=ounit, fmt=afmt, advance="no") v(i)
                enddo
                write(unit=ounit, fmt="(a)") ""
            endif
        endprocedure print_dvector

        module procedure print_ivector   
            integer                       :: m, i
            character(len=:), allocatable :: afmt
            logical                       :: aadvance

            if (present(advance)) then
                aadvance = advance
            else
                aadvance = .false.
            endif

            if (present(fmt)) then
                afmt = fmt
            else
                afmt = ivectorfmt
            endif

            m = size(v)

            if (present(message)) write (ounit, "(a)") message

            if (aadvance) then
                do i = 1, m
                    write(unit=ounit, fmt=afmt) v(i)
                enddo
            else
                do i = 1, m
                    write(unit=ounit, fmt=afmt, advance="no") v(i)
                enddo
                write(unit=ounit, fmt="(a)") ""
            endif
        endprocedure print_ivector
endsubmodule stduse_print