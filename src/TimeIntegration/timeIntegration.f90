#include "../macro.h"
module time_integration_m
    use sph
    use tools_m
    implicit none

    type Update
        real(8), allocatable :: Velocity(:)
        real(8) :: Density, Energy
    end type Update

    interface
        subroutine single_step(ntotal, ndummy, Particles, Delta, aver_v, Shear, dSdt)
            import Particle, Update
            integer,        intent(in)    :: ntotal
            integer,        intent(inout) :: ndummy
            type(Particle), intent(inout) :: Particles(:)
            type(Update),   intent(inout) :: Delta(:)
            real(8),        intent(inout) :: aver_v(:, :)
            real(8),        intent(in),    optional :: Shear(:, :, :)
            real(8),        intent(inout), optional :: dSdt(:, :, :)
        end subroutine single_step
    end interface

contains
    subroutine time_integration(ntotal, P)
        ! use ctrl_dict, only: Config
        integer, intent(in) :: ntotal
        type(Particle), intent(inout) :: P(:)
    
        call leap_frog(ntotal, P)
        
    end subroutine time_integration

    !!! Frog-Leap time integration
    subroutine leap_frog(ntotal, P)
        use ctrl_dict, only: Config, Field
        use cour_num_m
        use output_m
        implicit none
        integer, intent(in) :: ntotal
        type(Particle), intent(inout) :: P(:)
        integer :: ndummy = 0
        type(Update), allocatable :: Prev(:), Delta(:)
        real(8) :: aver_v(Field%dim, Field%maxn)
        real(8) :: time = 0
        real(8) :: aver_courant = 0, max_courant = 0, cntemp
#if SOLID
        real(8), dimension(Field%dim, Field%dim, Field%maxn) :: Shear_prev, Shear, dSdt
        real(8) :: J2, SigmaY
#endif

        integer i, k

        allocate(Prev(Field%maxn), Delta(Field%maxn))
        do i = 1, Field%maxn
            Prev(i)%Density = 0
            Prev(i)%Energy  = 0
            allocate(Prev(i)%Velocity(Field%dim), source=0._8)
            Delta(i)%Density = 0
            Delta(i)%Energy  = 0
            allocate(Delta(i)%Velocity(Field%dim), source=0._8)
        end do

        forall (i=1:Field%maxn)
            aver_v(:, i) = 0
#if SOLID
            Shear_prev(:, :, i) = 0
            Shear(:, :, i)      = 0
            dSdt(:, :, i)       = 0
#endif
        end forall

#if SOLID
        SigmaY = 6e8
#endif

#ifdef _OPENMP
        call omp_set_num_threads(Config%nthreads)
#endif

        call pbout(0, Config%max_time_step, .true.)
        call pbflush()

        do k = 1, Config%max_time_step
            Config%i_time_step = k

            if ( mod(Config%i_time_step, Config%print_interval) == 0 ) then
                call pbflush()
#ifndef _WIN32
                write(*, "(A)") repeat("—", 72)
#else
                write(*, "(A)") repeat("-", 72)
#endif
                write(*, "(2(A, G0))") " Courant Number mean: ", aver_courant, &
                                       " max: ", max_courant
                write(*,*) "Time step = ", to_string(Config%i_time_step)
                write(*, "(A, G0, A)") " deltaT = ", Config%delta_t, "s"
                write(*, "(A, G0, A)") " Time   = ", time + Config%delta_t, "s"
            end if

            !!! If not first time step, then update thermal energy, density
            !!! and velocity half a time step
            if ( Config%i_time_step /= 1 ) then
                !$OMP PARALLEL DO PRIVATE(i)
                do i = 1, ntotal
                    Prev(i)%Energy = P(i)%InternalEnergy
                    P(i)%InternalEnergy = P(i)%InternalEnergy + (Config%delta_t/2) * Delta(i)%Energy
                    if ( P(i)%InternalEnergy < 0 ) P(i)%InternalEnergy = 0
                    if ( .not. Config%sum_density_w ) then
                        Prev(i)%Density = P(i)%Density
                        P(i)%Density = P(i)%Density + (Config%delta_t/2) * Delta(i)%Density
                    end if
                    Prev(i)%Velocity = P(i)%v(:)
                    P(i)%v(:) = P(i)%v(:) + (Config%delta_t/2) * Delta(i)%Velocity(:)
#if SOLID
                    if ( abs(P(i)%Type) == 8 ) then
                        Shear_prev(:, :, i) = Shear(:, :, i)
                        Shear(:, :, i) = Shear(:, :, i) + (Config%delta_t/2)*dSdt(:, :, i)
                        J2 = sum( Shear(:, :, i)**2 )
                        Shear(:, :, i) = Shear(:, :, i) * min(1., sqrt(((SigmaY**2)/3)/J2))
                    end if
#endif
                end do
                !$OMP END PARALLEL DO
            end if

#if SOLID
            call single_step(ntotal, ndummy, P, Delta, aver_v, Shear, dSdt)
#else
            call single_step(ntotal, ndummy, P, Delta, aver_v)
#endif

            if ( Config%i_time_step == 1 ) then
                !$OMP PARALLEL DO PRIVATE(i)
                do i = 1, ntotal
                    P(i)%InternalEnergy = P(i)%InternalEnergy + (Config%delta_t/2) * Delta(i)%Energy
                    if ( P(i)%InternalEnergy < 0 ) P(i)%InternalEnergy = 0
                    if ( .not. Config%sum_density_w ) then
                        P(i)%Density = P(i)%Density + (Config%delta_t/2) * Delta(i)%Density
                    end if

                    P(i)%v(:) = P(i)%v(:) + (Config%delta_t/2) * Delta(i)%Velocity + aver_v(:, i)
                    P(i)%x(:) = P(i)%x(:) + Config%delta_t * P(i)%v(:)
#if SOLID
                    if ( abs(P(i)%Type) == 8 ) then
                        Shear(:, :, i) = Shear(:, :, i) + (Config%delta_t/2) * dSdt(:, :, i)
                        J2 = sum( Shear(:, :, i)**2 )
                        Shear(:, :, i) = Shear(:, :, i) * min(1., sqrt(((SigmaY**2)/3)/J2))
                    end if
#endif
                end do
                !$OMP END PARALLEL DO
            else
                max_courant = 0
                !$OMP PARALLEL DO PRIVATE(i) REDUCTION(max:max_courant) REDUCTION(+:aver_courant)
                do i = 1, ntotal
                    P(i)%InternalEnergy = Prev(i)%Energy + Config%delta_t * Delta(i)%Energy
                    if ( P(i)%InternalEnergy < 0 ) P(i)%InternalEnergy = 0
                    if ( .not. Config%sum_density_w ) then
                        P(i)%Density = Prev(i)%Density + Config%delta_t * Delta(i)%Density
                    end if
                    P(i)%v(:) = Prev(i)%Velocity + Config%delta_t * Delta(i)%Velocity + aver_v(:, i)
                    P(i)%x(:) = P(i)%x(:) + Config%delta_t * P(i)%v(:)
#if SOLID
                    if ( abs(P(i)%Type) == 8 ) then
                        Shear(:, :, i) = Shear_prev(:, :, i) + Config%delta_t * dSdt(:, :, i)
                        J2 = sum( Shear(:, :, i)**2 )
                        Shear(:, :, i) = Shear(:, :, i) * min(1., sqrt(((SigmaY**2)/3)/J2))
                    end if
#endif
                    cntemp = courant_num(P(i)%SmoothingLength, P(i)%divergenceVelocity, P(i)%SoundSpeed)
                    aver_courant = aver_courant + cntemp
                    if ( cntemp > max_courant ) max_courant = cntemp
                end do
                !$OMP END PARALLEL DO
                aver_courant = aver_courant / ntotal
            end if

            time = time + Config%delta_t

            if (mod(Config%i_time_step, Config%save_interval) == 0) then
                call output((Config%i_time_step/Config%save_interval), P(1:ntotal+ndummy))
            end if

            if ( mod(Config%i_time_step, Config%print_interval) == 0 ) then
                write(*,1000) "Location", "Velocity", "Acceleration"
                do i = 1, Field%dim
                    write(*,1001) P(Config%monitor_particle)%x(i), &
                                  P(Config%monitor_particle)%v(i), &
                                  Delta(Config%monitor_particle)%Velocity(i)
                end do
#ifndef _WIN32
                write(*, "(A)") repeat("—", 72)
#else
                write(*, "(A)") repeat("-", 72)
#endif
                write(*,*)
                call pbout(Config%i_time_step, Config%max_time_step, .true.)
            end if
        1000    format(1X, 3(A15))
        1001    format(1X, 3(2X, ES13.6))

        end do

    end subroutine leap_frog

end module time_integration_m