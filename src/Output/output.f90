module output_m
    use ctrl_dict, only: dim

    implicit none

contains
    subroutine output(index, ntotal, itype, x, v, mass, rho, p, e, c, hsml, div_r, Stress)
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
        real(8), intent(in) :: Stress(:, :, :)
        character(:), allocatable :: fileName
        ! logical :: types(8)
        integer :: particleType, num(-8:8)
        integer, allocatable :: type_list(:), type_indice(:, :)
        integer, dimension(ntotal) :: this_itype
        real(8), dimension(dim, ntotal) ::  this_x, this_v
        real(8), dimension(ntotal) :: this_mass, this_rho, this_p, this_e, &
                                    this_c, this_hsml, this_div_r
        real(8), dimension(dim, dim, ntotal) :: this_Stress
        

        integer i, j, k

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
                type_indice(5, num(5)+1:num(5)+num(0)) = type_indice(0, 1:num(0))
                num(5) = num(5) + num(0)
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
            do j = 1, num(particleType)
                k = type_indice(particleType, j)
                this_itype(j) = itype(k)
                this_x(:, j) = x(:, k)
                this_v(:, j) = v(:, k)
                this_mass(j) = mass(k)
                this_rho(j) = rho(k)
                this_p(j) = p(k)
                this_e(j) = e(k)
                this_c(j) = c(k)
                this_hsml(j) = hsml(k)
                this_div_r(j) = div_r(k)
                this_Stress(:, :, j) = Stress(:, :, k)
            end do
            call write_file(out_path//"/"//fileName//'.dat',               &
                            num(particleType), this_itype, this_x, this_v, &
                            this_mass, this_rho, this_p, this_e, this_c,   &
                            this_hsml, this_div_r, this_Stress)
            call write_vtk(vtk_path//"/"//fileName//'.vtk',               &
                           num(particleType), this_itype, this_x, this_v, &
                           this_mass, this_rho, this_p, this_e, this_c,   &
                           this_hsml, this_div_r, this_Stress)
            deallocate(fileName)
        end do

        ! call write_file(out_path//"/"//fileName//'.dat', &
        !                 ntotal, ndummy, itype, x, v, mass, rho, p, e, c, hsml, div_r)
        ! if ( write_vtk_w ) call write_vtk(vtk_path//"/"//"sph"//fileName//'.vtk', &
        !                                   ntotal, ndummy, itype, x, v, mass, rho, p, e, c, hsml, div_r)

    end subroutine output

    subroutine write_file(fileDir, ntotal, itype, x, v, mass, rho, p, e, c, hsml, div_r, Stress)
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
        real(8), intent(in) :: Stress(:, :, :)
        character(len=*), intent(in) :: fileDir
        integer i, d, dd

        open(11, file = fileDir)

        select case (dim)
        case (1)
            write(11, 1001) "Index", &
                            "X", "V", &
                            "Mass"           , "Density",     "Pressure", &
                            "InternalEnergy", "SoundSpeed", "Type",    "SmoothingLength", &
                            "DivDistance", &
                            "StressXX"
            do i = 1, ntotal
                write(11, 1002) i       , x(:, i), v(:, i), &
                                mass(i) , rho(i) , p(i)   , e(i), c(i), &
                                itype(i), hsml(i), div_r(i), ((Stress(dd, d, i), d = 1, dim), dd = 1, dim)
            end do
            1001 format(A5, 7(A17), A6, 3(A17))
            1002 format(I5, 7(2X, ES15.7E0), 2X, I4, 3(2X, ES15.7E0))

        case (2)
            write(11, 1003) "Index", &
                            "X", "Y", "U", "V", &
                            "Mass"           , "Density",     "Pressure", &
                            "InternalEnergy", "Sound Speed", "Type",    "SmoothingLength", &
                            "DivDistance", &
                            "StressXX", "StressXY", &
                            "StressYX", "StressYY"
            do i = 1, ntotal
                write(11, 1004) i       , x(:, i), v(:, i), &
                                mass(i) , rho(i) , p(i)   , e(i), c(i), &
                                itype(i), hsml(i), div_r(i), ((Stress(dd, d, i), d = 1, dim), dd = 1, dim)
            end do
            1003 format(A5, 9(A17), A6, 6(A17))
            1004 format(I5, 9(2X, ES15.7E0), 2X, I4, 6(2X, ES15.7E0))

        case (3)
            write(11, 1005) "Index", &
                            "X", "Y", "Z", "U", "V", "W", &
                            "Mass"           , "Density",     "Pressure", &
                            "InternalEnergy",  "SoundSpeed",  "Type",    "SmoothingLength", &
                            "DivDistance", &
                            "StressXX", "StressXY", "StressXZ", &
                            "StressYX", "StressYY", "StressYZ", &
                            "StressZX", "StressZY", "StressZZ"
            do i = 1, ntotal
                write(11, 1006) i       , x(:, i), v(:, i), &
                                mass(i) , rho(i) , p(i)   , e(i), c(i), &
                                itype(i), hsml(i), div_r(i), ((Stress(dd, d, i), d = 1, dim), dd = 1, dim)
            end do
            1005 format(A5, 11(A17), A6, 11(A17))
            1006 format(I5, 11(2X, ES15.7E0), 2X, I4, 11(2X, ES15.7E0))

        end select

        close(11)

    end subroutine write_file

    subroutine write_vtk(fileDir, ntotal, itype, x, v, mass, rho, p, e, c, hsml, div_r, Stress)
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
        real(8), intent(in) :: Stress(:, :, :)
        character(len=*), intent(in) :: fileDir
        integer i, d, dd

        open(11, file = fileDir)
        1001 format(*(ES12.5, 3X))

        !!! Write Header and particle coordinates
        write (11, "(A)") "# vtk DataFile Version 3.0, "//now()
        write (11, "(A)") "paraview_vtk_output"
        write (11, "(A)") "ASCII"
        write (11, "(A)") "DATASET UNSTRUCTURED_GRID"
        write (11, '(A, I0, A)') "POINTS ", ntotal, " double"
        do i = 1, ntotal
            write(11, 1001) x(:, i), (0.0, d = 1, 3-dim)
        end do
        write (11, "(A, I0)") "POINT_DATA ", ntotal

        !!! Write particle mass
        write (11, "(A)") "SCALARS Mass double 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) mass(i)
        end do

        !!! Write particle density
        write (11, "(A)") "SCALARS Density double 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) rho(i)
        end do

        !!! Write particle pressure
        write (11, "(A)") "SCALARS Pressure double 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) p(i)
        end do

        !!! Write particle Internal Energy
        write (11, "(A)") "SCALARS InternalEnergy double 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) e(i)
        end do

        !!! Write particle type
        write (11, "(A)") "SCALARS ParticleType int 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, "(I0)") itype(i)
        end do

        !!! Write particle smoothed length
        write (11, "(A)") "SCALARS SmoothingLength double 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) hsml(i)
        end do

        !!! Write particle distance divergence
        write (11, "(A)") "SCALARS DivDistance double 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) div_r(i)
        end do

        !!! Write particle velocity
        write (11, "(A)") "VECTORS U double"
        do i = 1, ntotal
            write(11, 1001) v(:, i), (0.0, d = 1, 3-dim)
        end do

        !!! Write particle stress
        write (11, "(A)") "TENSORS Stress double"
        do i = 1, ntotal
            write(11, 1001) ((Stress(d, dd, i), d = 1, dim), dd = 1, dim), (0.0, d = 1, 9-dim**2)
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

        1001 format(I6, 6(2X, ES15.7E0))
        1002 format(I6, 4(2X, ES15.7E0))
        1003 format(I6, 2X, I4, 2X, ES15.7E0)

        close(11)
        close(12)
        close(13)

    end subroutine output_3

end module output_m