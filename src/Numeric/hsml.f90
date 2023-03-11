module hsml_m
    use ctrl_dict, only: dim
    implicit none

contains
    !!! Subroutine to evolve smoothing length
    subroutine h_upgrade(ntotal, sle, delta_t, mass, rho, div_v, hsml)
        integer, intent(in) :: ntotal
        integer, intent(in) :: sle
        real(8), intent(in) :: delta_t
        real(8), intent(in) :: mass(:)
        real(8), intent(in) :: rho(:)
        real(8), intent(in) :: div_v(:)
        real(8), intent(inout) :: hsml(:)
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
            forall (i=1:ntotal) hsml(i) = factor * (mass(i)/rho(i))**(1._8/dim)
        case (2)
            !!! dh/dt = (-1/dim)*(h/rho)*(drho/dt)
            !!! drho/dt = sum(m*dv*dwdx )
            do i = 1, ntotal
                dhsmldt(i) = (hsml(i)/dim)*div_v(i)
                hsml(i)  = hsml(i) + delta_t*dhsmldt(i)
                if (hsml(i) <= 0) hsml(i) = hsml(i) - delta_t*dhsmldt(i)
            end do
        end select

    end subroutine h_upgrade

end module hsml_m