module hsml_m
    use ctrl_dict, only: dim
    use sph
    implicit none

contains
    !!! Subroutine to evolve smoothing length
    subroutine h_upgrade(ntotal, sle, delta_t, P)
        integer, intent(in) :: ntotal
        integer, intent(in) :: sle
        real(8), intent(in) :: delta_t
        type(Particle), intent(inout) :: P(:)
        real(8) :: factor
        real(8), allocatable :: dhsmldt(:)
        integer i

        allocate(dhsmldt(ntotal))

        select case (sle)
        case (0)
            !!! Keep smoothing length unchanged
            return
        case (1)
            factor = 2._8
            forall (i=1:ntotal) P(i)%SmoothingLength = factor*(P(i)%Mass/P(i)%Density)**(1._8/dim)
        case (2)
            !!! dh/dt = (-1/dim)*(h/rho)*(drho/dt)
            !!! drho/dt = sum(m*dv*dwdx )
            do i = 1, ntotal
                dhsmldt(i) = (P(i)%SmoothingLength/dim)*P(i)%divergenceVelocity
                P(i)%SmoothingLength  = P(i)%SmoothingLength + delta_t*dhsmldt(i)
                if (P(i)%SmoothingLength <= 0) then
                    P(i)%SmoothingLength = P(i)%SmoothingLength - delta_t*dhsmldt(i)
                end if
            end do
        end select

    end subroutine h_upgrade

end module hsml_m