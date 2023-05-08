module parse_toml_m
    use, intrinsic :: iso_fortran_env, only: err => error_unit
    use ctrl_dict, only: Config, Field, Project
    use tomlf, only: toml_table, get_value, toml_parse
    use tools_m, only: create_directory, ESC
    implicit none
    private

    public  :: fetch_control_value
contains
    subroutine fetch_control_value()

        ! call from_access_file()
        call get_project_dir_from_command_line()
        call get_details_from_sph_toml()

    end subroutine fetch_control_value

    subroutine get_project_dir_from_command_line()
        integer :: argument_num
        logical :: alive, flag = .false.
        character(512) :: buffer
        integer i

        argument_num = command_argument_count()
        select case( argument_num )
        case (0)
            write(*, "(A)") ESC//"[31m"
            write(*,*) "Please pass in the folder path as a command line argument."
            write(*, "(A)", advance="no") ESC//"[0m"
            stop
        case (1)
            call get_command_argument(1, buffer)
            allocate(Project%in_path, source=trim(adjustl(buffer)))
            inquire(file=Project%in_path//"/sph.toml", exist=alive)
            if ( .not. alive ) then
                write(*, "(A)") ESC//"[31m"
                write(*,*) "Target 'sph.toml' does not exist."
                write(*,*) "File path: "//Project%in_path//"/sph.toml"
                write(*, "(A)", advance="no") ESC//"[0m"
                stop
            end if
        case default
            do i = 1, argument_num
                call get_command_argument(i, buffer)
                allocate(Project%in_path, source=trim(adjustl(buffer)))
                inquire(file=Project%in_path//"/sph.toml", exist=alive)
                if ( alive ) then
                    flag = .true.
                    exit
                else
                    deallocate(Project%in_path)
                    cycle
                end if
            end do
            if ( flag ) then
                write(*, "(A)") ESC//"[33m"
                write(*,*) "Invalid argument ignored"
                write(*, "(A)", advance="no") ESC//"[0m"
            else
                write(*, "(A)") ESC//"[31m"
                write(*,*) "Could not find a valid path argument."
                write(*, "(A)", advance="no") ESC//"[0m"
                stop
            end if
        end select
        Project%out_path = Project%in_path//"/output"
        Project%vtk_path = Project%in_path//"/vtk"

        call create_directory(Project%out_path)
        call create_directory(Project%vtk_path)

    end subroutine get_project_dir_from_command_line

    subroutine get_details_from_sph_toml()
        !$ use omp_lib
        type(toml_table), allocatable :: sph_file
        type(toml_table), pointer :: subtable
        integer :: file_unit
#ifdef _OPENMP
        character(len=512) :: buffer
#endif

        open (newunit=file_unit, file=Project%in_path//"/sph.toml", status='old')
        call toml_parse(sph_file, file_unit)
        close (file_unit)

        call get_value(sph_file, 'name', Project%project_name, 'untitled')
        call get_value(sph_file, 'nick', Project%nick,         'unknown')
        write(*,*) ESC//"[32m"
#ifndef _WIN32
        write(*, "(A)") "╭"//repeat("─", 11 + len(Project%project_name))//"╮"
        write(*, "(A)") "│ Project: "//Project%project_name//" │"
        write(*, "(A)") "╰"//repeat("─", 11 + len(Project%project_name))//"╯"
#else
        write(*, "(A)") "|"//repeat("-", 11 + len(Project%project_name))//"|"
        write(*, "(A)") "| Project: "//Project%project_name//" |"
        write(*, "(A)") "|"//repeat("-", 11 + len(Project%project_name))//"|"
#endif
        write(*,*) ESC//"[0m"


        call get_value(sph_file, 'Parameter', subtable)

        call get_value(subtable, 'DIM',    Field%Dim)
        call get_value(subtable, 'maxN',   Field%Maxn,   5000)
        call get_value(subtable, 'ntotal', Field%ntotal      )
        call get_value(subtable, 'kPair',  Field%pairNum,  20  )

        call get_value(subtable, 'SPH',  Config%pa_sph, 2)
        call get_value(subtable, 'NNPS', Config%nnps,   1)
        call get_value(subtable, 'SKF',  Config%skf,    1)
        call get_value(subtable, 'SLE',  Config%sle,    0)

        write(*, "(A, I0)") " >> Nearest Neighbor Particle Search Method: ", Config%nnps
        write(*, "(A, I0)") " >> Smoothing Kernel Function: ",   Config%skf
        write(*, "(A, I0)") " >> Smoothing Length Estimation: ", Config%sle

        call get_value(subtable, 'monitorParticle', Config%monitor_particle, 1   )
        call get_value(subtable, 'printInterval',   Config%print_interval,   100 )
        call get_value(subtable, 'saveInterval',    Config%save_interval,    500 )
        call get_value(subtable, 'deltaT',          Config%delta_t,    dble(5e-4))
        call get_value(subtable, 'maxTimeStep',     Config%max_time_step,    1000)

        write(*,"(A, I0)") " >> Maximal Time Steps: ", Config%max_time_step


        call get_value(subtable, 'printStatistics',     Config%print_statistics_w, .true. )
        call get_value(subtable, 'summationDensity',    Config%sum_density_w,      .true. )
        call get_value(subtable, 'averageVelocity',     Config%aver_velocity_w,    .false.)
        call get_value(subtable, 'readInitialFile',     Config%read_ini_file_W,    .false.)
        call get_value(subtable, 'dummyParticle',       Config%dummy_parti_w,      .false.)
        call get_value(subtable, 'dpInput',             Config%dp_input,           .false.)
        call get_value(subtable, 'viscosity',           Config%viscosity_w,        .false.)
        call get_value(subtable, 'externalForce',       Config%ex_force_w,         .false.)
        call get_value(subtable, 'gravity',             Config%gravity_w,          .false.)
        call get_value(subtable, 'artificialViscosity', Config%arti_visc_w,        .false.)
        call get_value(subtable, 'artificialHeat',      Config%arti_heat_w,        .false.)
        call get_value(subtable, 'normalizeDensity',    Config%norm_dens_w,        .false.)
        call get_value(subtable, 'openBoundary',        Config%open_boundary_w,    .false.)

        if ( Config%dummy_parti_w .or. Config%gravity_w ) then
            Config%ex_force_w = .true.
        end if


        call get_value(sph_file, 'Parallel', subtable)

#ifdef _OPENMP
        Config%nthreads = 0

        call get_value(subtable, 'OMP_NUM_THREADS', Config%nthreads, 0)
        if ( Config%nthreads == 0 ) then 
            call get_environment_variable('OMP_NUM_THREADS', buffer)
            if ( trim(buffer) /= "" ) then
                read(buffer, *) Config%nthreads
            end if

            if ( Config%nthreads <= 0 ) then
                if ( Config%nthreads < 0) then
                    write(err, "(A, I0)") "Illegal number of threads: ", Config%nthreads
                end if
                Config%nthreads = omp_get_num_procs()
            end if
        end if

        write(*, "(A, I0)") " >> Number of Threads: ", Config%nthreads
#endif
        write(*,*)

        call get_value(sph_file, 'postProcess', subtable)

        call get_value(subtable, 'writeDummyParticle',     Config%write_dp_w,      .true.)
        call get_value(subtable, 'writeVTKfile',           Config%write_vtk_w,     .true.)
        call get_value(subtable, 'writeDummyParticle2VTK', Config%write_dp_vtk_w,  .true.)

        nullify (subtable)

    end subroutine get_details_from_sph_toml

    ! subroutine from_access_file()
    !     type(toml_table), allocatable :: access_file
    !     type(toml_table), pointer :: subtable
    !     integer file_unit

    !     open (newunit=file_unit, file='./access.toml', status='old')
    !     call toml_parse(access_file, file_unit)
    !     close (file_unit)

    !     call get_value(access_file, 'access', subtable)
    !     call get_value(subtable, 'Path',  Project%in_path,  './data/input')     !! Default Input  Dir
    !     Project%out_path = Project%in_path//"/output"
    !     Project%vtk_path = Project%in_path//"/vtk"

    !     nullify (subtable)

    ! end subroutine from_access_file

end module parse_toml_m
