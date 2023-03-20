module ctrl_dict
    ! use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
    ! use, intrinsic :: iso_c_binding, only: c_int, c_double
    implicit none
    private

    type Project_t
        character(:), allocatable, public :: in_path      !! Input  Directory
        character(:), allocatable, public :: out_path     !! Output Directory
        character(:), allocatable, public :: vtk_path     !! Output Directory
        character(:), allocatable, public :: project_name !! Current Project Name
        character(:), allocatable, public :: nick
    end type Project_t
    type(Project_t), public :: Project

    type Config_t
        !!! Time Integration Variables
        real(8) :: delta_t        !! Time step
        integer :: max_time_step  !! Maximum number of time steps to run
        integer :: i_time_step    !! Current time step
        !!! Output Parameters
        integer :: print_interval
        integer :: save_interval
        integer :: monitor_particle    !! Serial Number of the Monitored Particle
        !!! SPH Indicator Variety
        integer :: pa_sph  !! Indicator Variables for Particle Approximation
        integer :: nnps    !! Indicator Variables for the Nearest Neighbor Particle Search (NNPS)
        integer :: sle     !! Indicator Variables for Smooth Length Estimation
        integer :: skf     !! Indicator Variables for Smooth Kernel Function
        !!! Control Parameters
        logical :: sum_density_w    !! Whether to use Density Summation Method
        logical :: aver_velocity_w  !! Whether to use Speed Averaging Method
        logical :: read_ini_file_w  !! Control the Input of Initial Setup
        logical :: dummy_parti_w    !! Whether to Consider Dummy Particles
        logical :: dp_input         !! Control the Input of Dummy Particles
        logical :: viscosity_w      !! Whether to Consider Viscosity
        logical :: gravity_w        !! Whether to Consider Gravity
        logical :: ex_force_w       !! Whether to Consider External Force
        logical :: arti_visc_w      !! Whether to Consider Artificial Viscosity
        logical :: arti_heat_w      !! Whether to Consider Artificial Heat
        logical :: norm_dens_w      !! Whether to Normalize Density
        logical :: write_dp_w       !! Whether to write dummy particles
        logical :: write_vtk_w      !! Whether to write VTK files
        logical :: write_dp_vtk_w   !! Whether to write dummy particles to VTK files
        logical :: print_statistics_w  !! Control the Output of the SPH Particle Interaction State
        logical :: open_boundary_w

        !!! symmetry of the problem
        !!! nsym = 0 : no symmetry,
        !!!      = 1 : axis symmetry,
        !!!      = 2 : center symmetry.
        ! integer :: nsym
    
#ifdef _OPENMP
        integer :: nthreads
        integer :: chunkSize
#endif
    end type Config_t
    type(Config_t), public :: Config

    type Field_t
        !!! Physical Field Variables
        integer :: dim    !! Dimension
        integer :: maxn   !! Maximum Number of Particles
        integer :: ntotal
        ! real(8) :: tau_ba !! Viscous shear force

        !!! Particle Interaction Variables
        integer :: pairNum !! Maximum Allowed Number of Interacting Pairs

    end type Field_t
    type(Field_t), public :: Field

contains

end module ctrl_dict
