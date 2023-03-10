module tools_m
    use, intrinsic :: iso_fortran_env, only: stdout => output_unit
    implicit none
    character(len=*), parameter :: BS  = char(8)
    character(len=*), parameter :: CR  = char(13)
    character(len=*), parameter :: ESC = char(27)
    real(8), parameter :: PI = acos(-1.0_8)

    integer, private :: fdigit = 0

contains
    pure function to_string(var) result(str)
        implicit none
        class(*),     intent(in)  :: var
        character(:), allocatable :: str
        character(len=128)        :: buffer

        select type(var)
        type is (integer)
            write(buffer, *) var
            str = trim(adjustl(buffer))
        type is (real(8))
            write(buffer, *) var
            str = trim(adjustl(buffer))
        type is (real(4))
            write(buffer, *) var
            str = trim(adjustl(buffer))
        type is (logical)
            write(buffer, *) var
            str = trim(adjustl(buffer))
        end select

    end function to_string

    pure function round(var, decimal) result(r)
        implicit none
        class(*), intent(in) :: var
        integer,  intent(in) :: decimal
        integer :: aux, h
        real :: r

        select type(var)
        type is (real(8))
            associate (aux => mod(int(var*10.**(decimal+1)), 10))
                if ( aux < 5 ) then
                    r = floor(var*10.**(decimal+1) - aux) &
                        * 10.**(-decimal-1)
                else if ( aux == 5 ) then
                    associate (h => mod(int(var*10**(decimal+2)), 10))
                        if ( mod(h, 2) == 0 ) then
                            r = floor(var*10.**(decimal+1) - aux) &
                                * 10.**(-decimal-1)
                        else
                            r = floor(var*10.**(decimal+1) - aux + 10) &
                                * 10.**(-decimal-1)
                        end if
                    end associate
                else
                    r = floor(var*10.**(decimal+1) - aux + 10) &
                        * 10.**(-decimal-1)
                end if
            end associate
        type is (real(4))
            associate (aux => mod(int(var*10.**(decimal+1)), 10))
                if ( aux < 5 ) then
                    r = floor(var*10.**(decimal+1) - aux) &
                        * 10.**(-decimal-1)
                else if ( aux == 5 ) then
                    associate (h => mod(int(var*10.**(decimal+2)), 10))
                        if ( mod(h, 2) == 0 ) then
                            r = floor(var*10.**(decimal+1) - aux) &
                                * 10.**(-decimal-1)
                        else
                            r = floor(var*10.**(decimal+1) - aux + 10) &
                                * 10.**(-decimal-1)
                        end if
                    end associate
                else
                    r = floor(var*10.**(decimal+1) - aux + 10) &
                        * 10.**(-decimal-1)
                end if
            end associate
        end select

    end function round

    impure character(19) function now() result(t)
        character(len=8)  :: date
        character(len=10) :: time

        call date_and_time(date, time)

        t = date(1:4)//"-"// &
            date(5:6)//"-"// &
            date(7:8)//" "// &
            time(1:2)//":"// &
            time(3:4)//":"// &
            time(5:6)

    end function now

    impure subroutine print_error(var, info, type)
        implicit none
        class(*), intent(in) :: var
        character(len=*), intent(in) :: info, type

        write(*, "(A)", advance='no') ESC//"[31m"
        write(*, "(A)", advance='no') BS//type//" error:"
        write(*, "(A)", advance='no') ESC//"[0m"
        write(*,*) info//": "//to_string(var)//CR

    end subroutine print_error

    subroutine print_warning(var, info, type)
        implicit none
        class(*), intent(in) :: var
        character(len=*), intent(in) :: info, type

        write(*, "(A)", advance='no') ESC//"[33m"
        write(*, "(A)", advance='no') BS//type//" warning:"
        write(*, "(A)", advance='no') ESC//"[0m"
        write(*,*) info//": "//to_string(var)//CR

    end subroutine print_warning

    subroutine create_directory(name, info)
        implicit none
        character(len=*), intent(in) :: name
        character(len=*), intent(inout), optional :: info
        character(len=:), allocatable :: buffer
        logical :: exist

        allocate(buffer, source=trim(adjustl(name)))

#ifdef __INTEL_COMPILER
        inquire(Directory=buffer, Exist=exist)
#elif __GNUC__
        inquire(File=buffer, Exist=exist)
#endif
            if ( .not. exist ) then
                if ( present(info) ) then
                    info = "Created directory: "//trim(adjustl(buffer))
                end if
                call system('mkdir '//trim(adjustl(buffer)))
            else
                if ( present(info) ) then
                    info = "Directory already exists: "//trim(adjustl(buffer))
                end if
            end if
    
            deallocate(buffer)

    end subroutine create_directory

    !!! Thanks to Contrail (@StellaContrail) for subroutine pbout() and pbflush()
    !!! https://github.com/StellaContrail/FortranProgressbar.git
    !!! Modified by Teng Qing (@SZWHTQ)
    ! ????????????(Any)???????????????????????????????????????
    ! value / max : % (value out of max)
    ! isflushed : Flush????????????????????????????????????. ??????????????????????????????????????????????????????????????????????????????????????????. ?????????????????????false??????????????????
    subroutine pbout(value, max, isflushed)
        integer,intent(in) :: value, max
        double precision,save :: rate = 0d0, time = 0d0
        double precision dt, dr, estimate
        integer,parameter :: digit = 50
        integer remain(2), values(8), i
        ! character(len=20) :: FMT
        logical,intent(in), optional :: isflushed
        if (present(isflushed) .and. isflushed) then
            write (*, *)
        end if

        ! ??????max??????????????????????????????????????????????????????
        ! write (FMT, '(a, i0, a)') "(2(i", int(log10(real(max))) + 1, ",a))"

        ! ???????????????
        call date_and_time(values=values)
        dr = rate   ! ??????????????????rate???????????????dr???????????????
        dt = time   ! ??????????????????time???????????????dt???????????????

        ! ??????????????? : ?????????0??????????????????????????????????????????. -> ??????????????????????????????????????????????????????
        time = ((values(5)*60d0+values(6))*60d0+values(7))*1d3+values(8)
        dt = time - dt

        ! ???????????????
        rate = dble(value) / dble(max)
        dr = rate - dr

        ! ????????????????????? (milliseconds)
        estimate = (1d0 - rate) * (dt / dr)

        ! min/sec???????????????
        remain(1) = int(estimate/6d4) ! minutes
        remain(2) = int((estimate-remain(1)*6d4)*1d-3) ! seconds

        write (*, "(1X, 2(I0, A), $)") value, " / ", max, " ["
        do i = 1, int(digit*rate)
            write (*, '(a, $)') "="
        end do
        write (*, '(a, $)') ">"
        do i = int(digit*rate) + 1, digit
            write (*, '(a, $)') "-"
        end do

        if (isflushed) then
            fdigit = 2 * int(log10(real(max))+1) + 75
        end if
        write (*, '(a, f7.2, a, i3, a, i2, a, a, $)') &
               "]", 100d0*rate, "% ", remain(1), "m", remain(2), "s", CR

        if (value == max) then
            write (*, *)
        end if
    end subroutine

    ! ?????????????????????(Any)??????????????????????????????
    ! pbout()???isflushed???true???????????????????????????????????????false???????????????????????????
    subroutine pbflush()
        character(len=9) FMT
        if (fdigit == 0) then
            return
        end if

        write (FMT, '(a, i0, a)') "(", fdigit, "x, 3a)"
        write (*, FMT, advance='no') char(13), char(8), char(13)
    end subroutine

end module tools_m