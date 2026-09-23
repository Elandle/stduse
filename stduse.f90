module stduse
    use iso_fortran_env, only: real32, real64, input_unit, output_unit, iostat_end
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    implicit none

    integer, parameter :: sp     = real32
    integer, parameter :: dp     = real64
    integer, parameter :: stdin  = input_unit
    integer, parameter :: stdout = output_unit

    character(len=*), parameter :: dmatrixfmt = "(f17.8)"
    character(len=*), parameter :: dvectorfmt = "(f17.8)"
    character(len=*), parameter :: ivectorfmt = "(i6)"
    character(len=*), parameter :: imatrixfmt = "(i6)"

    interface print_matrix
        module procedure :: print_dmatrix
        module procedure :: print_imatrix
    endinterface

    interface print_vector
        module procedure :: print_dvector
        module procedure :: print_ivector
    endinterface

    abstract interface
        integer function seed_function(seed, i)
            integer, intent(in) :: seed
            integer, intent(in) :: i
        endfunction seed_function
    endinterface

    interface isclose
        module procedure :: isclose_dp
        module procedure :: isclose_zp
    endinterface isclose

    interface iscomplex
        module procedure :: iscomplex_zp
    endinterface iscomplex

    contains

        pure elemental integer function sgn(x)
            !
            ! Returns the sign of x as an integer.
            !
            !     1 if x >= 0
            !    -1 if x <  0
            !
            real(dp), intent(in) :: x

            if (x .ge. 0.0_dp) then
                sgn = 1
            else
                sgn = -1
            endif
        endfunction sgn

        pure elemental integer function del(i, j)
            integer, intent(in) :: i
            integer, intent(in) :: j

            if (i .eq. j) then
                del = 1
            else
                del = 0
            endif
        endfunction del

        !
        ! Sets seed to a "random" (at least according to Fortran's intrinsic
        ! random_init(repeatable=.false.)) single integer seed.
        !
        subroutine rand_seed(seed, sfunction)
            integer                 , intent(out) :: seed
            procedure(seed_function), optional    :: sfunction

            real(dp) :: r

            call random_init(repeatable=.false., image_distinct=.true.)
            call random_number(r)
            seed = int((2.0_dp*r - 1.0_dp) * huge(seed))
            call set_seed(seed, sfunction)
        end subroutine rand_seed

        !
        ! Initializes the random number generator using seed.
        ! There is an optional sfunction argument with signature:
        !
        !       sfunction(seed, i)
        !
        ! where seed and i are intent(in) integer and returns an integer.
        ! This function bridges the gap between single integer seeds
        ! and Fortran's seed array (it populates the seed array in a loop
        ! over this, ie Fortran's seed(i) = sfunction(seed, i))).
        ! A "random" one made by a bunch of bit operations is provided by default
        ! and does not overflow.
        !
        subroutine set_seed(seed, sfunction)
            integer                 , intent(in) :: seed
            procedure(seed_function), optional   :: sfunction

            integer              :: n
            integer              :: i
            integer, allocatable :: put(:)

            call random_seed(size=n)
            allocate(put(n))

            if (present(sfunction)) then
                do i = 1, n
                    put(i) = sfunction(seed, i)
                end do
            else
                do i = 1, n
                    put(i) = seed_function_dflt(seed, i)
                end do
            end if
            call random_seed(put=put)
            deallocate(put)
        endsubroutine set_seed

        function seed_function_dflt(seed, i) result(j)
            integer, intent(in) :: seed
            integer, intent(in) :: i

            integer :: j

            ! Bunch of random bit stuff (avoids overflow)
            j = ieor(seed, ishftc(i, 7))
            j = ieor(j, shiftl(j, 5))
            j = ieor(j, shiftr(j, 17))
            j = ishftc(j, 9)
            j = ieor(j, not(shiftr(j, 11)))
            j = ieor(j, shiftl(j, 3))
            j = ishftc(j, -7)
            j = ieor(j, ishftc(i, 13))
            j = ieor(j, shiftr(j, 19))
            j = ieor(j, shiftl(j, 6))
            j = ishftc(j, 11)
            j = ieor(j, shiftr(j, 8))
            j = ieor(j, shiftl(j, 15))
        endfunction seed_function_dflt

        subroutine print_dmatrix(A, ounit, message)
            real(dp)        , intent(in)           :: A(:, :)
            integer         , intent(in)           :: ounit
            character(len=*), intent(in), optional :: message
            
            integer :: m, n, i, j

            m = size(A, 1)
            n = size(A, 2)

            if (present(message)) then
                write (ounit, "(a)") message
            endif

            do i = 1, m
                do j = 1, n
                    write(unit=ounit, fmt=dmatrixfmt, advance="no") A(i, j)
                enddo
                write(unit=ounit, fmt="(a)") ""
            enddo
        endsubroutine print_dmatrix

        subroutine print_dvector(v, ounit, advance)
            real(dp), intent(in) :: v(:)
            integer , intent(in) :: ounit
            logical , optional   :: advance
            
            integer :: m, i
            logical :: adv

            if (present(advance)) then
                adv = advance
            else
                adv = .false.
            endif

            m = size(v)

            if (adv) then
                do i = 1, m
                    write(unit=ounit, fmt=dvectorfmt) v(i)
                enddo
            else
                do i = 1, m
                    write(unit=ounit, fmt=dvectorfmt, advance="no") v(i)
                enddo
                write(unit=ounit, fmt="(a)") ""
            endif
        endsubroutine print_dvector

        subroutine print_ivector(v, ounit, advance)
            integer, intent(in) :: v(:)
            integer, intent(in) :: ounit
            logical, optional   :: advance
            
            integer :: m, i
            logical :: adv

            if (present(advance)) then
                adv = advance
            else
                adv = .false.
            endif

            m = size(v)

            if (adv) then
                do i = 1, m
                    write(unit=ounit, fmt=ivectorfmt) v(i)
                enddo
            else
                do i = 1, m
                    write(unit=ounit, fmt=ivectorfmt, advance="no") v(i)
                enddo
                write(unit=ounit, fmt="(a)") ""
            endif
        endsubroutine print_ivector

        subroutine print_imatrix(A, ounit, message)
            integer         , intent(in)           :: A(:, :)
            integer         , intent(in)           :: ounit
            character(len=*), intent(in), optional :: message
            
            integer :: m, n, i, j

            m = size(A, 1)
            n = size(A, 2)

            if (present(message)) then
                write (ounit, "(a)") message
            endif

            do i = 1, m
                do j = 1, n
                    write(unit=ounit, fmt=imatrixfmt, advance="no") A(i, j)
                enddo
                write(ounit, "(a)") ""
            enddo
        endsubroutine print_imatrix

                !> \brief Reads the next line from a file.
        !!
        !! Reminder: in Fortran, files are streams.
        !! So if a line is read from a file (eg, by using the `read` intrinsic)
        !! that line is consumed in the stream: the next time a line is read,
        !! it is actually the next line in the file (to go back, an intrinsic
        !! such as `backspace` or `rewind` may be used).
        !!
        !! Reads the next line from file with unit `funit`, storing it in `line`
        !! with length `maxlen`. Leading spaces are removed, and `line` is padded
        !! with trailing spaces to fill up to `maxlen` length if necessary.
        !!
        !! Example: if the line being read is: </p>
        !! `nstab = 4` </p>
        !! Then if `readln` is called with `maxlen = 12`, the read in line is `"nstab = 4    "`. 
        !!
        !! \param[in]    funit  (`integer`)                      Unit of file to read from.
        !! \param[out]   iostat (`integer`)                      Status of the `read` statement on the file read from.
        !! \param[inout] line   (`character(len=:), allocatable`) String to allocate (`line` must be allocatable) read line into.
        !! \param[in]    maxlen (`integer, optional`)            Maximum length of a line to read and the exit length of `line`. Lines with data to be read in longer than `maxlen` being read into may result in errors. Default value: `100`.
        subroutine readln(funit, iostat, line, maxlen)
            integer                      , intent(in)              :: funit
            integer                      , intent(out)             :: iostat
            character(len=:), allocatable, intent(inout)           :: line
            integer                      , intent(in)   , optional :: maxlen

            ! Cannot edit maxlen since it is optional, so using another integer to hold the value used
            integer :: actual_maxlen

            if (present(maxlen)) then
                actual_maxlen = maxlen
            else
                ! default maximum length of 1024
                actual_maxlen = 1024
            endif

            ! allocate if not already allocated
            if (.not. allocated(line)) then
                allocate(character(actual_maxlen) :: line)
            endif

            ! force line to be maxlen long
            if (len(line) .ne. actual_maxlen) then
                deallocate(line)
                allocate(character(actual_maxlen) :: line)
            endif

            read(unit=funit, fmt="(a)", iostat=iostat) line
            line = ltrim(line)
        endsubroutine readln

        !> \brief Removes leading and trailing blank spaces from an input string.
        !!
        !! Example: if `line = "  abc de   "`, then `ltrim(line)` returns
        !! `"abc de"`.
        !!
        !! \param[in] line   (`character(len=*)`)       String to remove leading and trailing blank spaced from.
        !! \result    ltrim  (`character(len=trimlen)`) `line` with leading and trailing blank spaces removed, with a length `trimlen` of `line` minus the amount of leading and trailing blank spaces.
        function ltrim(line)
            character(len=*), intent(in) :: line

            character(len=len(trim(adjustl(line)))) :: ltrim

            ltrim = trim(adjustl(line))
        endfunction ltrim


        ! Fortran version of numpy's isclose
        ! https://numpy.org/doc/stable/reference/generated/numpy.isclose.html
        pure elemental function isclose_dp(a, b, rtol, atol, equal_nan) result(val)
            real(dp), intent(in)           :: a
            real(dp), intent(in)           :: b
            real(dp), intent(in), optional :: rtol
            real(dp), intent(in), optional :: atol
            logical,  intent(in), optional :: equal_nan

            logical  :: val
            real(dp) :: artol
            real(dp) :: aatol
            logical  :: aequal_nan
            logical  :: a_nan, b_nan

            ! Default values
            artol      = 1.0e-5_dp ! Relative tolerance
            aatol      = 1.0e-8_dp ! Absolute tolerance
            aequal_nan = .false.   ! Whether or not nan is treated as equal

            if (present(rtol))      artol      = rtol
            if (present(atol))      aatol      = atol
            if (present(equal_nan)) aequal_nan = equal_nan

            ! Check if input is nan
            a_nan = ieee_is_nan(a)
            b_nan = ieee_is_nan(b)
            if (a_nan .or. b_nan) then
                val = aequal_nan .and. a_nan .and. b_nan
                return
            endif

            ! Check actual equality
            if (a .eq. b) then
                val = .true.
                return
            endif

            ! Check if input is infinite
            if ((.not. ieee_is_finite(a)) .or. (.not. ieee_is_finite(b))) then
                val = .false.
                return
            end if

            ! Now guaranteed a and b are not nan or inf
            ! Check if they are within an absolute tolerance + relative tolerance.
            if (abs(a - b) .le. aatol + artol * abs(b)) then
                val = .true.
            else
                val = .false.
            endif
        end function isclose_dp

        pure elemental function isclose_zp(a, b, rtol, atol, equal_nan) result(val)
            complex(dp), intent(in)           :: a
            complex(dp), intent(in)           :: b
            real(dp)   , intent(in), optional :: rtol
            real(dp)   , intent(in), optional :: atol
            logical    , intent(in), optional :: equal_nan

            logical  :: val
            real(dp) :: artol
            real(dp) :: aatol
            logical  :: aequal_nan
            logical  :: a_nan, b_nan
            logical  :: a_finite, b_finite

            ! Default values
            artol      = 1.0e-5_dp ! Relative tolerance
            aatol      = 1.0e-8_dp ! Absolute tolerance
            aequal_nan = .false.   ! Whether or not nan is treated as equal

            if (present(rtol))      artol      = rtol
            if (present(atol))      aatol      = atol
            if (present(equal_nan)) aequal_nan = equal_nan

            ! Check if input is nan
            a_nan = ieee_is_nan(real(a, dp)) .or. ieee_is_nan(aimag(a))
            b_nan = ieee_is_nan(real(b, dp)) .or. ieee_is_nan(aimag(b))

            if (a_nan .or. b_nan) then
                val = aequal_nan .and. a_nan .and. b_nan
                return
            endif

            ! Check actual equality
            if (a .eq. b) then
                val = .true.
                return
            endif

            ! Check if input is finite
            a_finite = ieee_is_finite(real(a, dp)) .and. ieee_is_finite(aimag(a))
            b_finite = ieee_is_finite(real(b, dp)) .and. ieee_is_finite(aimag(b))

            if ((.not. a_finite) .or. (.not. b_finite)) then
                val = .false.
                return
            endif

            ! Now guaranteed a and b are not nan or inf
            ! Check if they are within an absolute tolerance + relative tolerance.
            if (abs(a - b) .le. aatol + artol * abs(b)) then
                val = .true.
            else
                val = .false.
            endif
        endfunction isclose_zp

        
        pure elemental function iscomplex_zp(x, rtol, atol) result(val)
            complex(dp), intent(in)           :: x
            real(dp)   , intent(in), optional :: rtol
            real(dp)   , intent(in), optional :: atol

            logical :: val

            val = .not. isclose(aimag(x), 0.0_dp, rtol=rtol, atol=atol)
        endfunction iscomplex_zp




endmodule stduse