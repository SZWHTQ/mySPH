program main
    use tools_m, only: print_error, print_warning

    implicit none
    integer :: dim = 4
    real(4) :: pi = 3.15
    real(8) :: dp = 3.142

    call print_error(dim, info="Dimension", type="value")
    call print_warning(pi, type="value",  info="Pi")
    call sub()
    call print_error(dp,  "Double precision", type="value")
    
end program main

subroutine sub()
    use tools_m, only: print_error

    implicit none
    real(8) :: dp = 3.142

    call print_error(dp, type="runtime", info="Double precision from sub")
    stop "at subroutine sub()"

end subroutine sub