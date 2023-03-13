module initial_m
    use ctrl_dict, only: maxn, max_interaction, dim
    use sph
    implicit none
    private
    
    integer, allocatable, public :: Type(:)  !! Type of Particles
    integer, allocatable, public :: status(:) !! Status of Particles
    real(8), allocatable, public :: x(:, :)   !! Position of Particles
    real(8), allocatable, public :: v(:, :)   !! Velocity of Particles
    real(8), allocatable, public :: mass(:)   !! Particle mass
    real(8), allocatable, public :: rho(:)    !! Particle density
    real(8), allocatable, public :: p(:)      !! Pressure
    real(8), allocatable, public :: e(:)      !! Internal energy
    real(8), allocatable, public :: c(:)      !! Speed of sound
    real(8), allocatable, public :: hsml(:)   !! Particle smooth length
    real(8), allocatable, public :: eta(:)    !! Shear viscosity coefficient
    integer, allocatable, public :: neighborList(:, :)   !! List of the Part of Interaction Pair
    real(8), allocatable, public :: w(:, :)         !! Smooth Kernel Function for a Given Interaction Pair
    real(8), allocatable, public :: dwdx(:, :, :)    !! The First Derivative of the Smooth Kernel Function for a Given Interaction Pair
    integer, allocatable, public :: neighborNum(:)

    public :: initialize

contains
    subroutine initialize()
        implicit none
        integer :: kpair

        kpair = max_interaction / maxn

        allocate(Type(maxn), status(maxn), source=0)
        allocate(x(dim, maxn), v(dim, maxn), source=0._8)
        allocate(mass(maxn), source=0._8)
        allocate(rho(maxn), p(maxn), e(maxn), c(maxn), source=0._8)
        allocate(hsml(maxn), eta(maxn), source=0._8)

        ! allocate(neighborList(max_interaction, 2), source=0)
        ! allocate(w(max_interaction), source=0._8)
        ! allocate(dwdx(dim, max_interaction), source=0._8)

        allocate(neighborNum(maxn), source=0)
        allocate(neighborList(maxn, kpair), source=0)
        allocate(w(maxn, kpair), source=0._8)
        allocate(dwdx(dim, maxn, kpair), source=0._8)

    end subroutine initialize


end module initial_m