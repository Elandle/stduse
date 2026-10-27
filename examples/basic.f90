program main
    use stduse
    implicit none

    real(dp) :: x(3)

    x(1) = 1.23_dp; x(2) = -2.561_dp; x(3) = 4.0_dp

    call print_vector(x, stdout)

endprogram main