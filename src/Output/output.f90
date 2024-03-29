module output_m
    use ctrl_dict, only: Project, Field

    implicit none

contains
    subroutine output(index, Particles, Delta)
        use sph
        ! use ctrl_dict,    only: write_vtk_w
        use tools_m,      only: to_string
        integer, intent(in) :: index
        type(Particle), intent(in) :: Particles(:)
        type(Update), intent(in), optional :: Delta(:)
        integer :: ntotal
        type(Particle), allocatable :: this(:)
        type(Update), allocatable :: thisDelta(:)
        integer :: particleType, num(-200:200)
        integer, allocatable :: type_list(:), type_indice(:, :)
        character(:), allocatable :: fileName

        integer i, j, k

        ntotal = size(Particles)
        call allocateParticleList(this, ntotal, Field%Dim, size(Particles(1)%neighborList))
        if ( present(Delta) ) then
            allocate(thisDelta(ntotal))
            do i = 1, ntotal
                thisDelta(i)%Density = 0
                thisDelta(i)%Energy  = 0
                allocate(thisDelta(i)%Velocity(Field%Dim), source=0._8)
            end do
        end if
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
                if ( present(Delta) ) then
                    thisDelta(j) = Delta(k)
                end if
            end do
            if ( present(Delta) ) then
                call write_file(Project%out_path//"/"//fileName//'.dat', this(1:num(particleType)), type_indice(particleType, :))
                call write_vtk(Project%vtk_path//"/"//fileName//'.vtk',  this(1:num(particleType)), thisDelta(1:num(particleType)))
            else
                call write_file(Project%out_path//"/"//fileName//'.dat', this(1:num(particleType)), type_indice(particleType, :))
                call write_vtk(Project%vtk_path//"/"//fileName//'.vtk',  this(1:num(particleType)))
            end if
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
        integer :: ntotal, fileUnit

        integer i

        ntotal = size(Particles)
        open(newunit=fileUnit, file = fileDir)

        if ( Particles(1)%Type <= 100 ) then
            select case (Field%Dim)
            case (1)
                write(fileUnit, 1001) "Index", "Type", "State", "X", "V",   &
                                "Mass", "Density",                          &
                                "Pressure", "InternalEnergy", "SoundSpeed", &
                                "SmoothingLength", "Viscosity",             &
                                "DivDistance",                              &
                                "Displacement_X",                           &
                                "StressXX"
                do i = 1, ntotal
                    if ( present(index) ) then
                        write(fileUnit, "(I8)", advance="no") index(i)
                    else
                        write(fileUnit, "(I8)", advance="no") i
                    end if
                    write(fileUnit, "(DT)") Particles(i)
                end do

            case (2)
                write(fileUnit, 1001) "Index", "Type", "State",             &
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
                        write(fileUnit, "(I8)", advance="no") index(i)
                    else
                        write(fileUnit, "(I8)", advance="no") i
                    end if
                    write(fileUnit, "(DT)") Particles(i)
                end do

            case (3)
                write(fileUnit, 1001) "Index", "Type", "State",             &
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
                        write(fileUnit, "(I8)", advance="no") index(i)
                    else
                        write(fileUnit, "(I8)", advance="no") i
                    end if
                    write(fileUnit, "(DT)") Particles(i)
                end do

            end select
        else
            select case (Field%Dim)
            case (1)
                write(fileUnit, 1001) "Index", "Type", "State", "X", "V",   &
                                "Mass", "Density",                          &
                                "Pressure", "InternalEnergy", "SoundSpeed", &
                                "SmoothingLength", "Viscosity",             &
                                "DivDistance",                              &
                                "Displacement_X",                           &
                                "StressXX", "vonMises"
                do i = 1, ntotal
                    if ( present(index) ) then
                        write(fileUnit, "(I8)", advance="no") index(i)
                    else
                        write(fileUnit, "(I8)", advance="no") i
                    end if
                    write(fileUnit, "(DT, ES24.16e3)") Particles(i), vonMises(Particles(i)%Stress)
                end do

            case (2)
                write(fileUnit, 1001) "Index", "Type", "State",             &
                                "X", "Y", "U", "V",                         &
                                "Mass" , "Density",                         &
                                "Pressure", "InternalEnergy", "SoundSpeed", &
                                "SmoothingLength", "Viscosity",             &
                                "DivDistance",                              &
                                "Displacement_X", "Displacement_Y",         &
                                "StressXX", "StressXY",                     &
                                "StressYX", "StressYY", "vonMises"
                do i = 1, ntotal
                    if ( present(index) ) then
                        write(fileUnit, "(I8)", advance="no") index(i)
                    else
                        write(fileUnit, "(I8)", advance="no") i
                    end if
                    write(fileUnit, "(DT, ES24.16e3)") Particles(i), vonMises(Particles(i)%Stress)
                end do

            case (3)
                write(fileUnit, 1001) "Index", "Type", "State",             &
                                "X", "Y", "Z", "U", "V", "W",               &
                                "Mass" , "Density",                         &
                                "Pressure", "InternalEnergy", "SoundSpeed", &
                                "SmoothingLength", "Viscosity",             &
                                "DivDistance",                              &
                                "Displacement_X", "Displacement_Y",         &
                                "Displacement_Z",                           &
                                "StressXX", "StressXY", "StressXZ",         &
                                "StressYX", "StressYY", "StressYZ",         &
                                "StressZX", "StressZY", "StressZZ", "vonMises"
                do i = 1, ntotal
                    if ( present(index) ) then
                        write(fileUnit, "(I8)", advance="no") index(i)
                    else
                        write(fileUnit, "(I8)", advance="no") i
                    end if
                    write(fileUnit, "(DT, ES24.16e3)") Particles(i), vonMises(Particles(i)%Stress)
                end do

            end select
        end if

        1001 format(3(A8), *(A17))

        close(fileUnit)

    end subroutine write_file

    subroutine write_vtk(fileDir, Particles, Delta)
        ! use ctrl_dict, only: write_dp_vtk_w
        use tools_m,   only: now
        use sph
        character(len=*), intent(in) :: fileDir
        type(Particle), intent(in) :: Particles(:)
        type(Update), intent(in), optional :: Delta(:)
        integer :: ntotal, fileUnit
        integer i, d, dd

        ntotal = size(Particles)
        open(newunit=fileUnit, file = fileDir)

        !!! Write Header and particle coordinates
        write (fileUnit, "(A)") "# vtk DataFile Version 3.0, "//now()
        write (fileUnit, "(A)") "paraview_vtk_output"
        write (fileUnit, "(A)") "ASCII"
        write (fileUnit, "(A)") "DATASET UNSTRUCTURED_GRID"
        write (fileUnit, '(A, I0, A)') "POINTS ", ntotal, " float"
        do i = 1, ntotal
            write(fileUnit, 1001) real(Particles(i)%x(:)), (0.0, d = 1, 3-Field%Dim)
        end do
        write (fileUnit, "(A, I0)") "POINT_DATA ", ntotal

        !!! Write particle mass
        write (fileUnit, "(A)") "SCALARS Mass float 1"
        write (fileUnit, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (fileUnit, 1001) real(Particles(i)%Mass)
        end do

        !!! Write particle density
        write (fileUnit, "(A)") "SCALARS Density float 1"
        write (fileUnit, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (fileUnit, 1001) real(Particles(i)%Density)
        end do

        !!! Write particle pressure
        write (fileUnit, "(A)") "SCALARS Pressure float 1"
        write (fileUnit, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (fileUnit, 1001) real(Particles(i)%Pressure)
        end do

        !!! Write particle Internal Energy
        write (fileUnit, "(A)") "SCALARS InternalEnergy float 1"
        write (fileUnit, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (fileUnit, 1001) real(Particles(i)%InternalEnergy)
        end do

        !!! Write particle type
        write (fileUnit, "(A)") "SCALARS Type int 1"
        write (fileUnit, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (fileUnit, "(I0)") Particles(i)%Type
        end do

        !!! Write particle state
        write (fileUnit, "(A)") "SCALARS State int 1"
        write (fileUnit, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (fileUnit, "(I0)") Particles(i)%State
        end do

        !!! Write particle state
        write (fileUnit, "(A)") "SCALARS Boundary int 1"
        write (fileUnit, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (fileUnit, "(I0)") Particles(i)%Boundary
        end do

        !!! Write particle smoothed length
        write (fileUnit, "(A)") "SCALARS SmoothingLength float 1"
        write (fileUnit, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (fileUnit, 1001) real(Particles(i)%SmoothingLength)
        end do

        !!! Write particle distance divergence
        write (fileUnit, "(A)") "SCALARS DivDistance float 1"
        write (fileUnit, "(A)") "LOOKUP_TABLE DEFAULT"
        do i = 1, ntotal
            write (fileUnit, 1001) real(Particles(i)%divergencePosition)
        end do

        !!! Write particle velocity
        write (fileUnit, "(A)") "VECTORS U float"
        do i = 1, ntotal
            write(fileUnit, 1001) real(Particles(i)%v(:)), (0.0, d = 1, 3-Field%Dim)
        end do

        !!! Write particle velocity
        write (fileUnit, "(A)") "VECTORS Displacement float"
        do i = 1, ntotal
            write(fileUnit, 1001) real(Particles(i)%Displacement(:)), (0.0, d = 1, 3-Field%Dim)
        end do

        !!! Write particle stress
        write (fileUnit, "(A)") "TENSORS Stress float"
        do i = 1, ntotal
            write(fileUnit, 1001) ((real(Particles(i)%Stress(d,dd)),d=1,Field%Dim),dd=1,Field%Dim), &
                            (0.0, d = 1, 9-Field%Dim**2)
        end do

        if ( Particles(1)%Type > 100 ) then
            !!! Write particle von Mises stress
            write (fileUnit, "(A)") "SCALARS vonMises float 1"
            write (fileUnit, "(A)") "LOOKUP_TABLE DEFAULT"
            do i = 1, ntotal
                write(fileUnit, 1001) real(vonMises(Particles(i)%Stress))
            end do
        end if

        if ( present(Delta) ) then
            !!! Write particle acceleration
            write (fileUnit, "(A)") "VECTORS Acceleration float"
            do i = 1, ntotal
                write(fileUnit, 1001) real(Delta(i)%Velocity(:)), (0.0, d = 1, 3-Field%Dim)
            end do
        end if

        1001 format(*(ES24.16e3, 3X))

        close(fileUnit)

    end subroutine write_vtk

    pure function vonMises(cauchyStress) result(Mises)
        real(8), intent(in) :: cauchyStress(:,:)
        real(8) :: Mises
        integer :: Dim

        Dim = size(cauchyStress,1)
        select case(Dim)
        case(1)
            Mises = abs(cauchyStress(1,1))
        case(2)
            Mises = sqrt((cauchyStress(1,1)+cauchyStress(2,2))**2 &
                            - 3*(cauchyStress(1,1)*cauchyStress(2,2)) &
                            + 3*(cauchyStress(1,2)**2))
        case(3)
            Mises = sqrt((cauchyStress(1,1)+cauchyStress(2,2)+cauchyStress(3,3))**2 &
                                - 3*(cauchyStress(2,2)*cauchyStress(3,3)  &
                                    +cauchyStress(3,3)*cauchyStress(1,1)  &
                                    +cauchyStress(1,1)*cauchyStress(2,2)) &
                                + 3*(cauchyStress(1,2)**2 &
                                    +cauchyStress(1,3)**2 &
                                    +cauchyStress(2,3)**2))
        end select

    end function vonMises

end module output_m