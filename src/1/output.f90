module output_m
    use ctrl_dict, only: dim

    implicit none

contains
    subroutine output(index, ntotal, itype, x, v, mass, rho, p, e, c, hsml, div_r)
        use parse_toml_m, only: out_path, vtk_path
        ! use ctrl_dict,    only: write_vtk_w
        use tools_m,      only: to_string
        integer, intent(in) :: index
        integer, intent(in) :: ntotal
        integer, intent(in) :: itype(:)
        real(8), intent(in) :: x(:, :)
        real(8), intent(in) :: v(:, :)
        real(8), intent(in) :: mass(:)
        real(8), intent(in) :: rho(:)
        real(8), intent(in) :: p(:)
        real(8), intent(in) :: e(:)
        real(8), intent(in) :: c(:)
        real(8), intent(in) :: hsml(:)
        real(8), intent(in) :: div_r(:)
        character(:), allocatable :: fileName
        ! logical :: types(8)
        integer :: particleType, num(-8:8)
        integer, allocatable :: type_list(:), type_indice(:, :)

        integer i

        ! types = .false.

        num = 0
        allocate(type_indice(-8:8, ntotal), source=0)
        do i = 1, ntotal
            particleType = itype(i)
            num(particleType) = num(particleType) + 1
            type_indice(particleType, num(particleType)) = i
        end do
        if ( num(0) /= 0 ) then
            if ( num(4) /= 0 ) then
                num(4) = num(4) + num(0)
                type_indice(4, num(4)+1:num(4)+num(0)) = type_indice(0, 1:num(0))
                type_indice(0, 1:num(0)) = 0
                num(0) = 0
            else
                num(5) = num(5) + num(0)
                type_indice(5, num(5)+1:num(5)+num(0)) = type_indice(0, 1:num(0))
                type_indice(0, 1:num(0)) = 0
                num(0) = 0
            end if
        end if

        allocate(type_list(0))
        do i = -8, 8
            if ( num(i) /= 0 ) type_list = [type_list, i]
        end do

        do i = 1, size(type_list)
            particleType = type_list(i)
            allocate(fileName, source="Type_"//to_string(particleType)//"_"//to_string(index))
            associate (slice => type_indice(particleType, 1:num(particleType)))
            call write_file(out_path//"/"//fileName//'.dat',                           &
                            num(particleType), itype(slice), x(:, slice), v(:, slice), &
                            mass(slice), rho(slice), p(slice), e(slice), c(slice), hsml(slice), div_r(slice))
            call write_vtk(vtk_path//"/"//fileName//'.vtk',                     &
                           num(particleType), itype(slice), x(:, slice), v(:, slice), &
                           mass(slice), rho(slice), p(slice), e(slice), c(slice), hsml(slice), div_r(slice))
            end associate
            deallocate(fileName)
        end do

        ! call write_file(out_path//"/"//fileName//'.dat', &
        !                 ntotal, ndummy, itype, x, v, mass, rho, p, e, c, hsml, div_r)
        ! if ( write_vtk_w ) call write_vtk(vtk_path//"/"//"sph"//fileName//'.vtk', &
        !                                   ntotal, ndummy, itype, x, v, mass, rho, p, e, c, hsml, div_r)

    end subroutine output

    subroutine write_file(fileDir, ntotal, itype, x, v, mass, rho, p, e, c, hsml, div_r)
        ! use ctrl_dict, only: write_dp_w
        integer, intent(in) :: ntotal
        integer, intent(in) :: itype(:)
        real(8), intent(in) :: x(:, :)
        real(8), intent(in) :: v(:, :)
        real(8), intent(in) :: mass(:)
        real(8), intent(in) :: rho(:)
        real(8), intent(in) :: p(:)
        real(8), intent(in) :: e(:)
        real(8), intent(in) :: c(:)
        real(8), intent(in) :: hsml(:)
        real(8), intent(in) :: div_r(:)
        character(len=*), intent(in) :: fileDir
        integer i

        open(11, file = fileDir)

        select case (dim)
        case (1)
            write(11, 1001) "Index", &
                            "X", "V", &
                            "Mass"           , "Density",     "Pressure", &
                            "Internal Energy", "Sound Speed", "Type",    "Smoothing Length", &
                            "DivDistance"
            do i = 1, ntotal
                write(11, 1002) i       , x(:, i), v(:, i), &
                                mass(i) , rho(i) , p(i)   , e(i), c(i), &
                                itype(i), hsml(i), div_r(i)
            end do
            1001 format(A5, 7(A17), A6, 2(A17))
            1002 format(I5, 7(2X, E15.8), 2X, I4, 2(2X, E15.8))

        case (2)
            write(11, 1003) "Index", &
                            "X", "Y", "U", "V", &
                            "Mass"           , "Density",     "Pressure", &
                            "Internal Energy", "Sound Speed", "Type",    "Smoothing Length", &
                            "DivDistance"
            do i = 1, ntotal
                write(11, 1004) i       , x(:, i), v(:, i), &
                                mass(i) , rho(i) , p(i)   , e(i), c(i), &
                                itype(i), hsml(i), div_r(i)
            end do
            1003 format(A5, 9(A17), A6, 2(A17))
            1004 format(I5, 9(2X, E15.8), 2X, I4, 2(2X, E15.8))

        case (3)
            write(11, 1005) "Index", &
                            "X", "Y", "Z", "U", "V", "W", &
                            "Mass"           , "Density",     "Pressure", &
                            "Internal Energy", "Sound Speed", "Type",    "Smoothing Length", &
                            "DivDistance"
            do i = 1, ntotal
                write(11, 1006) i       , x(:, i), v(:, i), &
                                mass(i) , rho(i) , p(i)   , e(i), c(i), &
                                itype(i), hsml(i), div_r(i)
            end do
            1005 format(A5, 11(A17), A6, 2(A17))
            1006 format(I5, 11(2X, E15.8), 2X, I4, 2(2X, E15.8))

        end select

        close(11)

    end subroutine write_file

    subroutine write_vtk(fileDir, ntotal, itype, x, v, mass, rho, p, e, c, hsml, div_r)
        ! use ctrl_dict, only: write_dp_vtk_w
        use tools_m,   only: now
        integer, intent(in) :: ntotal
        integer, intent(in) :: itype(:)
        real(8), intent(in) :: x(:, :)
        real(8), intent(in) :: v(:, :)
        real(8), intent(in) :: mass(:)
        real(8), intent(in) :: rho(:)
        real(8), intent(in) :: p(:)
        real(8), intent(in) :: e(:)
        real(8), intent(in) :: c(:)
        real(8), intent(in) :: hsml(:)
        real(8), intent(in) :: div_r(:)
        character(len=*), intent(in) :: fileDir
        integer i

        open(11, file = fileDir)
        1001 format(*(ES12.5, 3X))

        !!! Write Header and particle coordinates
        write (11, "(A)") "# vtk DataFile Version 3.0, "//now()
        write (11, "(A)") "paraview_vtk_output"
        write (11, "(A)") "ASCII"
        write (11, "(A)") "DATASET UNSTRUCTURED_GRID"
        write (11, '(A, I0, A)') "POINTS ", ntotal, " float"
        do i = 1, ntotal
            select case(dim)
            case (1)
                write(11, 1001) x(:, i), 0.0, 0.0
            case (2)
                write(11, 1001) x(:, i), 0.0
            case (3)
                write(11, 1001) x(:, i)
            end select
        end do
        write (11, "(A, I0)") "POINT_DATA ", ntotal

        !!! Write particle mass
        write (11, "(A)") "SCALARS Mass float 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) mass(i)
        end do

        !!! Write particle density
        write (11, "(A)") "SCALARS Density float 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) rho(i)
        end do

        !!! Write particle pressure
        write (11, "(A)") "SCALARS Pressure float 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) p(i)
        end do

        !!! Write particle Internal Energy
        write (11, "(A)") "SCALARS Internal_Energy float 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) e(i)
        end do

        !!! Write particle type
        write (11, "(A)") "SCALARS Particle_Type int 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, "(I0)") itype(i)
        end do

        !!! Write particle smoothed length
        write (11, "(A)") "SCALARS Smoothing_Length float 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) hsml(i)
        end do

        !!! Write particle distance divergence
        write (11, "(A)") "SCALARS DivDistance float 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) div_r(i)
        end do

        !!! Write particle velocity
        write (11, "(A)") "VECTORS U float"
        do i = 1, ntotal
            select case(dim)
            case (1)
                write(11, 1001) v(:, i), 0.0, 0.0
            case (2)
                write(11, 1001) v(:, i), 0.0
            case (3)
                write(11, 1001) v(:, i)
            end select
        end do

    end subroutine write_vtk


    subroutine output_3(ntotal, ndummy, itype, x, v, mass, rho, p, e, hsml)
        use ctrl_dict,    only: i_time_step
        use parse_toml_m, only: out_path
        use tools_m,      only: to_string
        integer, intent(in) :: ntotal
        integer, intent(in) :: ndummy
        integer, intent(in) :: itype(:)
        real(8), intent(in) :: x(:, :)
        real(8), intent(in) :: v(:, :)
        real(8), intent(in) :: mass(:)
        real(8), intent(in) :: rho(:)
        real(8), intent(in) :: p(:)
        real(8), intent(in) :: e(:)
        real(8), intent(in) :: hsml(:)
        integer i

        open(11, file = out_path // '/f_'//to_string(i_time_step)//'xv.dat')
        open(12, file = out_path // '/f_'//to_string(i_time_step)//'state.dat')
        open(13, file = out_path // '/f_'//to_string(i_time_step)//'other.dat')

        write(11, *) to_string(ntotal+ndummy)
        do i = 1, ntotal + ndummy
            write(11, 1001) i, x(:, i),  v(:, i)
            write(12, 1002) i, mass(i),  rho(i), p(i), e(i)
            write(13, 1003) i, itype(i), hsml(i)
        end do

        1001 format(I6, 6(2X, E15.8))
        1002 format(I6, 4(2X, E15.8))
        1003 format(I6, 2X, I4, 2X, E15.8)

        close(11)
        close(12)
        close(13)

    end subroutine output_3

end module output_m