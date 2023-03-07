program main
    use tools_m, only: to_string
    implicit none
    real(8) :: PI_dp = acos(-1._8)
    real(4) :: PI_sp = acos(-1._4)
    integer :: PI_int = 3
    logical :: flag = .true. 

    write(*,*) "start "//to_string(PI_dp)//" end"
    write(*,*) "start "//to_string(PI_sp)//" end"
    write(*,*) "start "//to_string(PI_int)//" end"
    write(*,*) "start "//to_string(flag)//" end"
    
end program main