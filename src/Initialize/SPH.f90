#include "../macro.h"
module sph
    implicit none
    private

    type, public :: Particle
        integer :: Type, State
        real(8), allocatable :: x(:), v(:)
        real(8) :: Mass, Density
        real(8) :: Pressure, InternalEnergy, SoundSpeed, SmoothingLength, Viscocity
        real(8) :: divergencePosition, divergenceVelocity
        real(8), allocatable :: Stress(:, :)
        integer :: neighborNum
        integer, allocatable :: neighborList(:)
        real(8), allocatable :: w(:)
        real(8), allocatable :: dwdx(:, :)
    contains
        private
        procedure, pass :: allocate => allocateParticle
        procedure, pass :: write => writeParticle
        procedure, pass :: read  => readParticle
        generic, public :: write(formatted) => write
        generic, public :: read(formatted)  => read
    end type Particle

    public :: allocateParticleList, allocateNeighborList

contains
    pure subroutine allocateParticle(P, Dim, PairNum)
        class(Particle), intent(inout) :: P
        integer, intent(in) :: Dim, PairNum

        allocate(P%x(Dim), P%v(Dim), source=0.0_8)
        allocate(P%Stress(Dim, Dim), source=0.0_8)
        allocate(P%neighborList(PairNum), source=0)
        allocate(P%w(PairNum), source=0.0_8)
        allocate(P%dwdx(Dim, PairNum), source=0.0_8)

    end subroutine allocateParticle

    subroutine writeParticle(P, unit, iotype, v_list, iostat, iomsg)
        class(Particle), intent(in) :: P
        integer, intent(in) :: unit
        character(len=*), intent(in) :: iotype
        integer, intent(in), optional :: v_list(:)
        integer, intent(out), optional :: iostat
        character(len=*), intent(inout), optional :: iomsg

        integer, parameter :: scalarNum = 8, vectorNum = 2, tensorNum = 1
        character(len=128) :: fmt
        integer i, j, dim, N

        dim = size(P%x)

        N = scalarNum + vectorNum * dim + tensorNum * dim * dim

#if WRITE_NEIGHBOR_INFO
        write(fmt, "(A, 3(I0, A))") "(2(I8), 2X", N, "(ES15.7e0, 2X), "  &
                                // "I0, 2X, ", P%neighborNum, "(I0, 2X)", &
                                P%neighborNum + P%neighborNum * dim, "(ES15.7e0, 2X))"
        write(unit, fmt) P%Type, P%State,                           &
                         P%x, P%v, P%Mass, P%Density,                &
                         P%Pressure, P%InternalEnergy, P%SoundSpeed, &
                         P%SmoothingLength, P%Viscocity,      &
                         P%divergencePosition,                       &
                         ((P%Stress(i, j), j=1, dim), i=1, dim),     &
                         P%neighborNum, P%neighborList, P%w, P%dwdx
#else
        write(fmt, "(A, 3(I0, A))") "(2(I8), 2X, ", N, "(ES15.7e0, 2X))"
        write(unit, fmt) P%Type, P%State,                          &
                         P%x, P%v, P%Mass, P%Density,                &
                         P%Pressure, P%InternalEnergy, P%SoundSpeed, &
                         P%SmoothingLength, P%Viscocity,      &
                         P%divergencePosition,                       &
                         ((P%Stress(i, j), j=1, dim), i=1, dim)
#endif
    end subroutine writeParticle

    subroutine readParticle(P, unit, iotype, v_list, iostat, iomsg)
        class(Particle), intent(inout) :: P
        integer, intent(in) :: unit
        character(len=*), intent(in), optional :: iotype
        integer, intent(in), optional :: v_list(:)
        integer, intent(out), optional :: iostat
        character(len=*), intent(inout), optional :: iomsg

        integer i, j, dim

        dim = size(P%x)

#if WRITE_NEIGHBOR_INFO
        read(unit, *) P%Type, P%State,                          &
                      P%x, P%v, P%Mass, P%Density,                &
                      P%Pressure, P%InternalEnergy, P%SoundSpeed, &
                      P%SmoothingLength, P%Viscocity,      &
                      P%divergencePosition,                       &
                      ((P%Stress(i, j), j=1, dim), i=1, dim),     &
                      P%neighborNum, P%neighborList, P%w, P%dwdx
#else
        read(unit, *) P%Type, P%State,                          &
                      P%x, P%v, P%Mass, P%Density,                &
                      P%Pressure, P%InternalEnergy, P%SoundSpeed, &
                      P%SmoothingLength, P%Viscocity,      &
                      P%divergencePosition,                       &
                      ((P%Stress(i, j), j=1, dim), i=1, dim)
#endif
    end subroutine readParticle

    subroutine allocateParticleList(Particles, ParticleNum, Dim, PairNum)
        type(Particle), allocatable, intent(inout) :: Particles(:)
        integer, intent(in) :: ParticleNum, Dim, PairNum

        integer i

        allocate(Particles(ParticleNum))

        !$omp parallel do private(i)
        do i = 1, ParticleNum
            Particles(i)%Type               = 0
            Particles(i)%State              = 0
            Particles(i)%Mass               = 0.0_8
            Particles(i)%Density            = 0.0_8
            Particles(i)%Pressure           = 0.0_8
            Particles(i)%InternalEnergy     = 0.0_8
            Particles(i)%SoundSpeed         = 0.0_8
            Particles(i)%SmoothingLength    = 0.0_8
            Particles(i)%Viscocity   = 0.0_8
            Particles(i)%divergencePosition = 0.0_8
            Particles(i)%divergenceVelocity = 0.0_8
            allocate(Particles(i)%x(Dim), Particles(i)%v(Dim), source=0.0_8)
            allocate(Particles(i)%Stress(Dim, Dim), source=0.0_8)
        end do
        !$omp end parallel do

        call allocateNeighborList(Particles, Dim, PairNum)

    end subroutine allocateParticleList

    subroutine allocateNeighborList(Particles, Dim, PairNum)
        type(Particle), intent(inout) :: Particles(:)
        integer, intent(in) :: Dim, PairNum

        integer i

        !$omp parallel do private(i)
        do i = 1, size(Particles)
            if ( Particles(i)%State /= 0 ) cycle
            Particles(i)%neighborNum = 0
            if ( allocated(Particles(i)%neighborList) ) then
                deallocate(Particles(i)%neighborList)
            end if
            allocate(Particles(i)%neighborList(PairNum), source=0)
            if ( allocated(Particles(i)%w) ) then
                deallocate(Particles(i)%w)
            end if
            allocate(Particles(i)%w(PairNum), source=0.0_8)
            if ( allocated(Particles(i)%dwdx) ) then
                deallocate(Particles(i)%dwdx)
            end if
            allocate(Particles(i)%dwdx(Dim, PairNum), source=0.0_8)
        end do
        !$omp end parallel do

    end subroutine allocateNeighborList

end module sph