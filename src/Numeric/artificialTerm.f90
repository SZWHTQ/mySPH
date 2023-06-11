module AT_m !!! Artificial Term
    use SPH
    use ctrl_dict, only: Config, Field, Project
    implicit none

contains
    subroutine artificialStress(ntotal, P, dvdt)
        use kernel_m, only: kernel
        use tools_m, only: PI
        integer, intent(in) :: ntotal
        type(Particle), intent(inout) :: P(:)
        real(8), intent(inout) :: dvdt(:, :)
        ! integer :: N
        real(8) :: theta, c, s, aux
        real(8), allocatable :: principalStress(:), R(:, :, :), Temp(:)
        real(8) :: delta, wd, F
        real(8), parameter :: lambda = 0.2, exponential = 4

        integer i, j, k, d, e

        ! N = size(P)

        do i = 1, ntotal
            dvdt(:, i) = 0
        end do

        select case(Project%nick)
        case ("waterImpact")
            delta = 0.002
        case ("can_beam")
            delta = 1e-3
        end select

        allocate(principalStress(Field%Dim), R(Field%Dim, Field%Dim, ntotal), source=0._8)
        allocate(Temp(Field%Dim), source=0._8)

        select case (Field%Dim)
        case (1)
            error stop "Artificial stress for 1D problem is still incomplete."
        case (2)
            !$omp parallel do private(i, j, k, d, theta, c, s, aux, principalStress, Temp)
            do i = 1, ntotal
                ! if ( P(i)%Type < 100 ) then
                !     cycle
                ! end if
                aux = (P(i)%Stress(1,1)-P(i)%Stress(2,2))
                if (aux == 0 ) then
                    theta = 0.5 * PI
                else
                    theta = 0.5 * atan((P(i)%Stress(1,2)+P(i)%Stress(2,1))/aux)
                end if
                s = sin(theta)
                c = cos(theta)
                principalStress(1) = c**2 * P(i)%Stress(1,1) + s**2 * P(i)%Stress(2,2)  &
                                       + 2 * s * c * P(i)%Stress(1,2)
                principalStress(2) = s**2 * P(i)%Stress(1,1) + c**2 * P(i)%Stress(2,2)  &
                                       - 2 * s * c * P(i)%Stress(1,2)
                do d = 1, Field%Dim
                    if ( principalStress(d) > 0 ) then
                        Temp(d) = -lambda * principalStress(d) / (P(i)%Density**2)
                    else
                        Temp(d) = 0
                    end if
                end do
                R(1, 1, i) = c**2 * Temp(1) + s**2 * Temp(2)
                R(2, 2, i) = s**2 * Temp(1) + c**2 * Temp(2)
                R(1, 2, i) = s * c * (Temp(1) - Temp(2))
                R(2, 1, i) = R(1, 2, i)
            end do
            !$omp end parallel do

            !$omp parallel do private(i, j, k, d, e, wd, F)
            do i = 1, ntotal
                do k = 1, P(i)%neighborNum
                    j = P(i)%neighborList(k)
                    ! if ( P(j)%Type < 100 ) then
                    !     cycle
                    ! end if

                    Temp = 0
                    call kernel(delta, 1*Temp, P(i)%SmoothingLength, wd, Temp)
                    F = (P(i)%w(k) / wd) ** exponential

                    do d = 1, Field%Dim
                        do e = 1, Field%Dim
                            dvdt(d, i) = dvdt(d, i) + &
                                P(j)%Mass * (R(d, e, i) + R(d, e, j)) * F * P(i)%dwdx(e, k)
                        end do
                    end do

                end do
            end do
            !$omp end parallel do
        case (3)
            error stop "Artificial stress for 3D problem is still incomplete."
        end select

    end subroutine artificialStress

end module AT_m