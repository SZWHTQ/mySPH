program main
    use parse_toml_m
    use initial_m
    use input_m
    use time_integration_m

    implicit none
    integer startT, endT, rate
    character(len=16) :: buffer

    call system_clock(startT)

    call fetch_control_value()

    call initialize()

    call input()

    call time_integration()

    call system_clock(endT, rate)
    write(buffer, "(F15.3)") dble((endT - startT))/rate
    write(*,*) "Time used: ", trim(adjustl(buffer)), "s"


end program main
