program main
    use tools_m, only: round
    implicit none
    real(8) :: PI_8 = 3.151
    real(4) :: PI_4 = 3.152


    write(*, "(F8.1)") PI_8
    write(*, "(F8.1)") PI_4

    write(*, "(F8.1)") round(PI_8, 1)
    write(*, "(F8.1)") round(PI_4, 1)

    write(*, "(F8.4)") round(3.14159, 4)

end program main