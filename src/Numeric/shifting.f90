module shifting_m
    use SPH
    use ctrl_dict, only: Config, Field
    implicit none

contains
    subroutine shifting(ntotal, P)
        integer, intent(in) :: ntotal
        type(Particle), intent(inout) :: P(:)
        real(8), allocatable :: concentration(:), shift(:)
        real(8) :: diffusion
        real(8), parameter :: lambda = 0.1

        integer i, j, k

        allocate(concentration(Field%Dim), shift(Field%Dim), source=0._8)
        
        do i = 1, ntotal
            if ( p(i)%Boundary == 1 ) then
                cycle
            end if
            concentration = 0
            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)
                concentration(:) = concentration(:) + P(j)%Mass/P(j)%Density * P(i)%dwdx(:, k)
            end do
            diffusion = lambda * P(i)%SmoothingLength**2! / Config%delta_t
            shift = -diffusion * concentration! * Config%delta_t
            P(i)%x = P(i)%x + shift
        end do

    end subroutine shifting

end module shifting_m