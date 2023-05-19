module output_m
    use ctrl_dict, only: Project, Field

    implicit none

contains
    subroutine output(index, Particles)
        use sph
        ! use ctrl_dict,    only: write_vtk_w
        use tools_m,      only: to_string
        integer, intent(in) :: index
        type(Particle), intent(inout) :: Particles(:)
        integer :: ntotal
        type(Particle), allocatable :: this(:)
        integer :: particleType, num(-200:200)
        integer, allocatable :: type_list(:), type_indice(:, :)
        character(:), allocatable :: fileName

        integer i, j, k

        ntotal = size(Particles)
        call allocateParticleList(this, ntotal, Field%Dim, size(Particles(1)%neighborList))
        num = 0
        allocate(type_indice(-200:200, ntotal), source=0)
        do i = 1, ntotal
            particleType = Particles(i)%Type
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
        do i = -200, 200
            if ( num(i) /= 0 ) type_list = [type_list, i]
        end do

        do i = 1, size(type_list)
            particleType = type_list(i)
            allocate(fileName, source="Type_"//to_string(particleType)//"_"//to_string(index))
            do j = 1, num(particleType)
                k = type_indice(particleType, j)
                this(j) = Particles(k)
            end do
            call write_file(Project%out_path//"/"//fileName//'.dat', this(1:num(particleType)), type_indice(particleType, :))
            call write_vtk(Project%vtk_path//"/"//fileName//'.vtk',  this(1:num(particleType)))
            deallocate(fileName)
        end do

        ! call write_file(Project%out_path//"/"//to_string(index)//'.dat', Particles)

        deallocate(this, type_list, type_indice)

    end subroutine output

    subroutine write_file(fileDir, Particles, index)
        use sph
        character(len=*), intent(in) :: fileDir
        type(Particle), intent(in) :: Particles(:)
        integer, intent(in), optional :: index(:)
        integer :: ntotal

        integer i

        ntotal = size(Particles)
        open(11, file = fileDir)

        select case (Field%Dim)
        case (1)
            write(11, 1001) "Index", "Type", "State", "X", "V",         &
                            "Mass", "Density",                          &
                            "Pressure", "InternalEnergy", "SoundSpeed", &
                            "SmoothingLength", "Viscosity",             &
                            "DivDistance",                              &
                            "Displacement_X",                           &
                            "StressXX"
            do i = 1, ntotal
                if ( present(index) ) then
                    write(11, "(I8)", advance="no") index(i)
                else
                    write(11, "(I8)", advance="no") i
                end if
                write(11, "(DT)") Particles(i)
            end do
            1001 format(3(A8), 12(A17))

        case (2)
            write(11, 1002) "Index", "Type", "State",                   &
                            "X", "Y", "U", "V",                         &
                            "Mass" , "Density",                         &
                            "Pressure", "InternalEnergy", "SoundSpeed", &
                            "SmoothingLength", "Viscosity",             &
                            "DivDistance",                              &
                            "Displacement_X", "Displacement_Y",         &
                            "StressXX", "StressXY",                     &
                            "StressYX", "StressYY"
            do i = 1, ntotal
                if ( present(index) ) then
                    write(11, "(I8)", advance="no") index(i)
                else
                    write(11, "(I8)", advance="no") i
                end if
                write(11, "(DT)") Particles(i)
            end do
            1002 format(3(A8), 18(A17))

        case (3)
            write(11, 1003) "Index", "Type", "State",                   &
                            "X", "Y", "Z", "U", "V", "W",               &
                            "Mass" , "Density",                         &
                            "Pressure", "InternalEnergy", "SoundSpeed", &
                            "SmoothingLength", "Viscosity",             &
                            "DivDistance",                              &
                            "Displacement_X", "Displacement_Y",         &
                            "Displacement_Z",                           &
                            "StressXX", "StressXY", "StressXZ",         &
                            "StressYX", "StressYY", "StressYZ",         &
                            "StressZX", "StressZY", "StressZZ"
            do i = 1, ntotal
                if ( present(index) ) then
                    write(11, "(I8)", advance="no") index(i)
                else
                    write(11, "(I8)", advance="no") i
                end if
                write(11, "(DT)") Particles(i)
            end do
            1003 format(3(A8), 25(A17))

        end select

        close(11)

    end subroutine write_file

    subroutine write_vtk(fileDir, Particles)
        ! use ctrl_dict, only: write_dp_vtk_w
        use tools_m,   only: now
        use sph
        character(len=*), intent(in) :: fileDir
        type(Particle), intent(in) :: Particles(:)
        integer :: ntotal
        integer i, d, dd

        ntotal = size(Particles)
        open(11, file = fileDir)
        1001 format(*(ES24.16e3, 3X))

        !!! Write Header and particle coordinates
        write (11, "(A)") "# vtk DataFile Version 3.0, "//now()
        write (11, "(A)") "paraview_vtk_output"
        write (11, "(A)") "ASCII"
        write (11, "(A)") "DATASET UNSTRUCTURED_GRID"
        write (11, '(A, I0, A)') "POINTS ", ntotal, " float"
        do i = 1, ntotal
            write(11, 1001) real(Particles(i)%x(:)), (0.0, d = 1, 3-Field%Dim)
        end do
        write (11, "(A, I0)") "POINT_DATA ", ntotal

        !!! Write particle mass
        write (11, "(A)") "SCALARS Mass float 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) real(Particles(i)%Mass)
        end do

        !!! Write particle density
        write (11, "(A)") "SCALARS Density float 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) real(Particles(i)%Density)
        end do

        !!! Write particle pressure
        write (11, "(A)") "SCALARS Pressure float 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) real(Particles(i)%Pressure)
        end do

        !!! Write particle Internal Energy
        write (11, "(A)") "SCALARS InternalEnergy float 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) real(Particles(i)%InternalEnergy)
        end do

        !!! Write particle type
        write (11, "(A)") "SCALARS Type int 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, "(I0)") Particles(i)%Type
        end do

        !!! Write particle state
        write (11, "(A)") "SCALARS State int 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, "(I0)") Particles(i)%State
        end do

        !!! Write particle smoothed length
        write (11, "(A)") "SCALARS SmoothingLength float 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) real(Particles(i)%SmoothingLength)
        end do

        !!! Write particle distance divergence
        write (11, "(A)") "SCALARS DivDistance float 1"
        write (11, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (11, 1001) real(Particles(i)%divergencePosition)
        end do

        !!! Write particle velocity
        write (11, "(A)") "VECTORS U float"
        do i = 1, ntotal
            write(11, 1001) real(Particles(i)%v(:)), (0.0, d = 1, 3-Field%Dim)
        end do

        !!! Write particle velocity
        write (11, "(A)") "VECTORS Displacement float"
        do i = 1, ntotal
            write(11, 1001) real(Particles(i)%Displacement(:)), (0.0, d = 1, 3-Field%Dim)
        end do

        !!! Write particle stress
        write (11, "(A)") "TENSORS Stress float"
        do i = 1, ntotal
            write(11, 1001) ((real(Particles(i)%Stress(d,dd)),d=1,Field%Dim),dd=1,Field%Dim), &
                            (0.0, d = 1, 9-Field%Dim**2)
        end do

    end subroutine write_vtk

end module output_m