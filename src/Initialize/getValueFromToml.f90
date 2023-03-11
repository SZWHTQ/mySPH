module parse_toml_m
    use ctrl_dict
    use tomlf, only: toml_table, get_value, toml_parse
    implicit none
    private
    character(len=*), parameter :: ESC = char(27)

    character(:), allocatable, public :: in_path      !! Input  Directory
    character(:), allocatable, public :: out_path     !! Output Directory
    character(:), allocatable, public :: vtk_path     !! Output Directory
    character(:), allocatable, public :: project_name !! Current Project Name
    character(:), allocatable, public :: nick

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

        argument_num = iargc()
        select case( argument_num )
        case (0)
            write(*, "(A)") ESC//"[31m"
            write(*,*) "Please pass in the folder path as a command line argument."
            write(*, "(A)", advance="no") ESC//"[0m"
            stop
        case (1)
            call getarg(1, buffer)
            allocate(in_path, source=trim(adjustl(buffer)))
            inquire(file=in_path//"/sph.toml", exist=alive)
            if ( .not. alive ) then
                write(*, "(A)") ESC//"[31m"
                write(*,*) "Target 'sph.toml' does not exist."
                write(*,*) "File path: "//in_path//"/sph.toml"
                write(*, "(A)", advance="no") ESC//"[0m"
                stop
            end if
        case default
            do i = 1, argument_num
                call getarg(i, buffer)
                allocate(in_path, source=trim(adjustl(buffer)))
                inquire(file=in_path//"/sph.toml", exist=alive)
                if ( alive ) then
                    flag = .true.
                    exit
                else
                    deallocate(in_path)
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
        out_path = in_path//"/output"
        vtk_path = in_path//"/vtk"


    end subroutine get_project_dir_from_command_line

    subroutine get_details_from_sph_toml()
        !$ use omp_lib
        type(toml_table), allocatable :: sph_file
        type(toml_table), pointer :: subtable
        integer :: kpair
        integer :: file_unit
#ifdef _OPENMP
        character(len=512) :: buffer
#endif

        open (newunit=file_unit, file=in_path//"/sph.toml", status='old')
        call toml_parse(sph_file, file_unit)
        close (file_unit)

        call get_value(sph_file, 'name', project_name, 'untitled')
        call get_value(sph_file, 'nick', nick,         'unknown')
        write(*,*) ESC//"[32m"
        write(*, "(A)") "╭"//repeat("─", 11 + len(project_name))//"╮"
        write(*, "(A)") "│ Project: "//project_name//" │"
        write(*, "(A)") "╰"//repeat("─", 11 + len(project_name))//"╯"
        write(*,*) ESC//"[0m"


        call get_value(sph_file, 'Parameter', subtable)

        call get_value(subtable, 'DIM',    dim)
        call get_value(subtable, 'maxN',   maxn,   5000)
        call get_value(subtable, 'ntotal', ntotal      )
        call get_value(subtable, 'kPair',  kpair,  20  )
        max_interaction = kpair * maxn

        call get_value(subtable, 'SPH',  pa_sph, 2)
        call get_value(subtable, 'NNPS', nnps,   1)
        call get_value(subtable, 'SKF',  skf,    1)
        call get_value(subtable, 'SLE',  sle,    0)

        write(*, "(A, I0)") " >> Nearest Neighbor Particle Search Method: ", nnps
        write(*, "(A, I0)") " >> Smoothing Kernel Function: ",   skf
        write(*, "(A, I0)") " >> Smoothing Length Estimation: ", sle

        call get_value(subtable, 'monitorParticle', monitor_particle, 1   )
        call get_value(subtable, 'printInterval',   print_interval,   100 )
        call get_value(subtable, 'saveInterval',    save_interval,    500 )
        call get_value(subtable, 'deltaT',          delta_t,    dble(5e-4))
        call get_value(subtable, 'maxTimeStep',     max_time_step,    1000)

        write(*,"(A, I0)") " >> Maximal Time Steps: ", max_time_step


        call get_value(subtable, 'printStatistics',     print_statistics_w, .true. )
        call get_value(subtable, 'summationDensity',    sum_density_w,      .true. )
        call get_value(subtable, 'averageVelocity',     aver_velocity_w,    .false.)
        call get_value(subtable, 'readInitialFile',     read_ini_file_W,    .false.)
        call get_value(subtable, 'dummyParticle',       dummy_parti_w,      .false.)
        call get_value(subtable, 'dpInput',             dp_input,           .false.)
        call get_value(subtable, 'viscosity',           viscosity_w,        .false.)
        call get_value(subtable, 'externalForce',       ex_force_w,         .false.)
        call get_value(subtable, 'gravity',             gravity_w,          .false.)
        call get_value(subtable, 'artificialViscosity', arti_visc_w,        .false.)
        call get_value(subtable, 'artificialHeat',      arti_heat_w,        .false.)
        call get_value(subtable, 'normalizeDensity',    norm_dens_w,        .false.)
        call get_value(subtable, 'DSPH',                DSPH_w,             .false.)

        call get_value(subtable, 'symmetry', nsym, 0)

        call get_value(subtable, 'waterEOS', eos_water_form, 2)

        if ( dummy_parti_w .or. gravity_w ) then
            ex_force_w = .true.
        end if


        call get_value(sph_file, 'Parallel', subtable)

#ifdef _OPENMP
        call getenv('OMP_NUM_THREADS', buffer)
        if ( trim(buffer)/="" ) then
            read(buffer, *) nthreads
            call get_value(subtable, 'OMP_NUM_THREADS', nthreads, 1*nthreads)
        else
            call get_value(subtable, 'OMP_NUM_THREADS', nthreads, 1)
        end if
        if ( nthreads == 0 ) then
            if ( trim(buffer)/="" ) then
                read(buffer, *) nthreads
            else
                nthreads = omp_get_num_procs()
            end if
        end if
        write(*, "(A, I0)") " >> Number of Threads: ", nthreads
#endif
        write(*,*)

        call get_value(sph_file, 'postProcess', subtable)

        call get_value(subtable, 'writeDummyParticle',     write_dp_w,      .true.)
        call get_value(subtable, 'writeVTKfile',           write_vtk_w,     .true.)
        call get_value(subtable, 'writeDummyParticle2VTK', write_dp_vtk_w,  .true.)

        nullify (subtable)

    end subroutine get_details_from_sph_toml

    subroutine from_access_file()
        type(toml_table), allocatable :: access_file
        type(toml_table), pointer :: subtable
        integer file_unit

        open (newunit=file_unit, file='./access.toml', status='old')
        call toml_parse(access_file, file_unit)
        close (file_unit)

        call get_value(access_file, 'access', subtable)
        call get_value(subtable, 'Path',  in_path,  './data/input')     !! Default Input  Dir
        out_path = in_path//"/output"
        vtk_path = in_path//"/vtk"

        nullify (subtable)

    end subroutine from_access_file

end module parse_toml_m
