submodule(stduse) stduse_utilities
    implicit none

    contains

        module procedure sgn_sp
            sgn_sp = int(sign(1.0_sp, x))
        endprocedure sgn_sp

        module procedure sgn_dp
            sgn_dp = int(sign(1.0_dp, x))
        endprocedure sgn_dp

        module procedure sgn_i
            sgn_i = sign(1, x)
        endprocedure sgn_i

        module procedure del
            del = merge(1, 0, i .eq. j)
        endprocedure del

        ! Fortran version of numpy's isclose
        ! https://numpy.org/doc/stable/reference/generated/numpy.isclose.html
        module procedure isclose_dp
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
        endprocedure isclose_dp

        module procedure isclose_zp
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
        endprocedure isclose_zp

        module procedure iscomplex_zp
            val = .not. isclose(aimag(x), 0.0_dp, rtol=rtol, atol=atol)
        endprocedure iscomplex_zp
endsubmodule stduse_utilities