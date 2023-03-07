program main
    use, intrinsic :: iso_fortran_env, only: rk => real32, stdout => output_unit
    use tools_m, only: pbflush, pbout, to_string
    implicit none
    integer :: i, N = 5, sleepTime = 1

    do i = 1, N
        call pbflush()
        ! write(*,*) i
        ! write(*, '(a,i0)') 'Current number of time step = ', i
        ! write(*, "(a,g0.3,a)") 'Total time = ', N, 's'
        write(*,*) "DEMO: ", to_string(i)
        call pbout(i, N, .true.)
        call sleep(sleepTime)
    end do

    stop
end program main