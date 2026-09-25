module stduse
    use iso_fortran_env, only: real32, real64, input_unit, output_unit, iostat_end
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
    implicit none
    private
    ! ||||||||||||||||||||||||||||||
    ! ||| Publicity declarations |||
    ! ||||||||||||||||||||||||||||||

    ! constants -----------------------------------------------------------------------------------
    public :: dp
    public :: sp
    public :: stdin
    public :: stdout
    ! ---------------------------------------------------------------------------------------------

    ! character -----------------------------------------------------------------------------------
    public :: ltrim
    public :: readln
    ! ---------------------------------------------------------------------------------------------

    ! print ---------------------------------------------------------------------------------------
    public :: print_matrix
    public :: print_vector
    ! ---------------------------------------------------------------------------------------------

    ! random --------------------------------------------------------------------------------------
    public  :: rand_seed
    public  :: seed_function
    private :: seed_function_dflt
    public  :: set_seed
    ! ---------------------------------------------------------------------------------------------

    ! utilities -----------------------------------------------------------------------------------
    public :: del
    public :: isclose
    public :: iscomplex
    public :: sgn
    ! ---------------------------------------------------------------------------------------------


    ! ||||||||||||||||||||||||||||||||||||||||||
    ! ||| Interface and constant definitions |||
    ! ||||||||||||||||||||||||||||||||||||||||||

    ! constants -----------------------------------------------------------------------------------
    integer, parameter :: sp     = real32
    integer, parameter :: dp     = real64
    integer, parameter :: stdin  = input_unit
    integer, parameter :: stdout = output_unit
    ! ---------------------------------------------------------------------------------------------

    ! character -----------------------------------------------------------------------------------
    interface ltrim
        module function ltrim(line)
            character(len=*), intent(in) :: line
            character(len=:), allocatable :: ltrim
        endfunction ltrim
    endinterface ltrim
    interface readln
        module subroutine readln(funit, iostat, line, maxlen)
            integer                      , intent(in)              :: funit
            integer                      , intent(out)             :: iostat
            character(len=:), allocatable, intent(inout)           :: line
            integer                      , intent(in)   , optional :: maxlen
        endsubroutine readln
    endinterface readln
    ! ---------------------------------------------------------------------------------------------

    ! print ---------------------------------------------------------------------------------------
    interface print_matrix
        module subroutine print_dmatrix(A, ounit, message, fmt)
            real(dp)        , intent(in)           :: A(:, :)
            integer         , intent(in)           :: ounit
            character(len=*), intent(in), optional :: message
            character(len=*), intent(in), optional :: fmt
        endsubroutine print_dmatrix
        module subroutine print_imatrix(A, ounit, message, fmt)
            integer         , intent(in)           :: A(:, :)
            integer         , intent(in)           :: ounit
            character(len=*), intent(in), optional :: message
            character(len=*), intent(in), optional :: fmt
        endsubroutine print_imatrix
    endinterface print_matrix
    interface print_vector
        module subroutine print_dvector(v, ounit, message, advance, fmt)
            real(dp)        , intent(in)           :: v(:)
            integer         , intent(in)           :: ounit
            character(len=*), intent(in), optional :: message
            logical         , intent(in), optional :: advance
            character(len=*), intent(in), optional :: fmt
        endsubroutine print_dvector
        module subroutine print_ivector(v, ounit, message, advance, fmt)
            integer         , intent(in)           :: v(:)
            integer         , intent(in)           :: ounit
            character(len=*), intent(in), optional :: message
            logical         , intent(in), optional :: advance
            character(len=*), intent(in), optional :: fmt
        endsubroutine print_ivector
    endinterface print_vector
    ! ---------------------------------------------------------------------------------------------

    ! random --------------------------------------------------------------------------------------
    interface rand_seed
        module subroutine rand_seed(seed, sfunction)
            integer                 , intent(out) :: seed
            procedure(seed_function), optional    :: sfunction
        endsubroutine rand_seed
    endinterface rand_seed
    abstract interface
        integer function seed_function(seed, i)
            integer, intent(in) :: seed
            integer, intent(in) :: i
        endfunction seed_function
    endinterface
    interface seed_function_dflt
        module function seed_function_dflt(seed, i) result(j)
            integer, intent(in) :: seed
            integer, intent(in) :: i
            integer :: j
        endfunction seed_function_dflt
    endinterface seed_function_dflt
    interface set_seed
        module subroutine set_seed(seed, sfunction)
            integer                 , intent(in) :: seed
            procedure(seed_function), optional   :: sfunction
        endsubroutine set_seed
    endinterface set_seed
    ! ---------------------------------------------------------------------------------------------

    ! utilities -----------------------------------------------------------------------------------
    interface del
        pure elemental integer module function del(i, j)
            integer, intent(in) :: i
            integer, intent(in) :: j
        endfunction del
    endinterface del
    interface isclose
        pure elemental module function isclose_dp(a, b, rtol, atol, equal_nan) result(val)
            real(dp), intent(in)           :: a
            real(dp), intent(in)           :: b
            real(dp), intent(in), optional :: rtol
            real(dp), intent(in), optional :: atol
            logical,  intent(in), optional :: equal_nan
            logical  :: val
        endfunction isclose_dp
        pure elemental module function isclose_zp(a, b, rtol, atol, equal_nan) result(val)
            complex(dp), intent(in)           :: a
            complex(dp), intent(in)           :: b
            real(dp)   , intent(in), optional :: rtol
            real(dp)   , intent(in), optional :: atol
            logical    , intent(in), optional :: equal_nan
            logical  :: val
        endfunction isclose_zp
    endinterface isclose
    interface iscomplex
        pure elemental module function iscomplex_zp(x, rtol, atol) result(val)
            complex(dp), intent(in)           :: x
            real(dp)   , intent(in), optional :: rtol
            real(dp)   , intent(in), optional :: atol
            logical :: val
        endfunction iscomplex_zp
    endinterface iscomplex
    interface sgn
        pure elemental integer module function sgn_dp(x)
            real(dp), intent(in) :: x
        endfunction sgn_dp
        pure elemental integer module function sgn_i(x)
            integer, intent(in) :: x
        endfunction sgn_i
        pure elemental integer module function sgn_sp(x)
            real(sp), intent(in) :: x
        endfunction sgn_sp
    endinterface sgn
    ! ---------------------------------------------------------------------------------------------
endmodule stduse