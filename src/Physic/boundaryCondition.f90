module boundary_condition_m
    use ctrl_dict, only: Project, Config, Field
    use sph,       only: Particle
    use geometry_m
    implicit none

    private

contains
    subroutine gen_non_reflecting_bc(ntotal, Particles)
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: Particles(:)

        select case (Project%nick)
        case ("undex_chamber")
            call gen_undex_chamber_bc(ntotal, Particles)
        end select

    end subroutine gen_non_reflecting_bc

    subroutine calc_nrbc_property(Particles)
        type(Particle), intent(inout) :: Particles(:)

        select case (Project%nick)
        case ("undex_chamber")
            call calc_undex_chamber_nrbc_property(Particles)
        end select

    end subroutine calc_nrbc_property

    !!! Open Boundary Non-Reflecting Condition
    subroutine gen_undex_chamber_nrbc(ntotal, Particles)
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: Particles(:)
        integer :: nbuffer, nghost
        type(rectangle_t) :: FluidDomain, AllDomain
        type(point_t) :: point
        real(8), dimension(2, 4) :: fixedPoints, normal
        integer :: entry_count = 0
        real(8) :: dx
        save FluidDomain, fixedPoints, entry_count, dx

        integer :: N, layer = 4
        integer i, k, l

        nbuffer = 0
        nghost = 0
        entry_count = entry_count + 1
        if ( entry_count == 1 ) then
            dx = abs(Particles(2)%x(1) - Particles(1)%x(1))
            FluidDomain = rectangle_t(   &
                [real(8) :: -0.5, -0.5], &
                [real(8) ::  1.0,  1.0], &
                0 )
            AllDomain = rectangle_t(     &
                [real(8) :: -0.5, -0.5], &
                ([real(8) :: 1.0,  1.0]  &
                    + 2 * layer * dx)  , &
                0 )

            !! Left
            fixedPoints(:, 1) = FluidDomain%center        &
                - [real(8) :: FluidDomain%length(1)/2, 0] &
                - [real(8) :: dx, 0] / 2
            normal(:, 1) = [1, 0]

            !! Bottom
            fixedPoints(:, 2) = FluidDomain%center        &
                - [real(8) :: 0, FluidDomain%length(2)/2] &
                - [real(8) :: 0, dx] / 2
            normal(:, 2) = [0, 1]

            !! Right
            fixedPoints(:, 3) = FluidDomain%center        &
                + [real(8) :: FluidDomain%length(1)/2, 0] &
                + [real(8) :: dx, 0] / 2
            normal(:, 3) = [-1, 0]

            !! Top
            fixedPoints(:, 4) = FluidDomain%center        &
                + [real(8) :: 0, FluidDomain%length(2)/2] &
                + [real(8) :: 0, dx] / 2
            normal(:, 4) = [0, -1]

            !!! This N is for particles' number in each layer
            N = int(1 / dx) + 1
            do l = 1, layer
                !! Left
                do i = 1, N
                    nbuffer = nbuffer + 1
                    Particles(ntotal+nbuffer)%x = [-0.5,  0.5] + [- l, 1 - i] * dx
                    Particles(ntotal+nbuffer)%State = 1
                end do
                !! Bottom
                do i = 1, N
                    nbuffer = nbuffer + 1
                    Particles(ntotal+nbuffer)%x = [-0.5, -0.5] + [1 + i, - l] * dx
                    Particles(ntotal+nbuffer)%State = 1
                end do
                !! Right
                do i = 1, N
                    nbuffer = nbuffer + 1
                    Particles(ntotal+nbuffer)%x = [ 0.5, -0.5] + [  l, 1 + i] * dx
                    Particles(ntotal+nbuffer)%State = 1
                end do
                !! Top
                do i = 1, N
                    nbuffer = nbuffer + 1
                    Particles(ntotal+nbuffer)%x = [ 0.5,  0.5] + [1 - i,   l] * dx
                    Particles(ntotal+nbuffer)%State = 1
                end do
            end do
            nghost = nbuffer
            do i = ntotal + nbuffer + 1, ntotal + nbuffer + nghost
                Particles(i)%x &
                    = calcGhostPosition(point_t(Particles(i - nghost)%x, 0))
                Particles(i)%State = 2
                Particles(i - nghost)%neighborNum = 1
                Particles(i - nghost)%neighborList(1) = i
            end do
        else
            do i = 1, ntotal
                point = point_t(Particles(i)%x, 0)
                if ( Particles(i)%State == 0 ) then !! Fluid Particle
                    if ( BufferDomainContain(point) ) then
                        !! Fluid Particle converts to Buffer Particle
                        Particles(i)%State = 1
                        nghost = nghost + 1
                        Particles(ntotal+nbuffer+nghost)%x &
                            = calcGhostPosition(point_t(Particles(i)%x, 0))
                        Particles(ntotal+nbuffer+nghost)%State = 2
                        Particles(i)%neighborNum = 1
                        Particles(i)%neighborList(1) = ntotal+nbuffer+nghost
                    end if !! BufferDomainContain(point)
                else if ( Particles(i)%State == 1 ) then !! Buffer Particle
                    if ( .not. BufferDomainContain(point) ) then
                        if ( FluidDomain%contain(point) ) then
                            !! Buffer Particle converts to Fluid Particle
                            Particles(i)%State = 0
                            Particles(Particles(i)%neighborList(1))%State = -1
                            !! Add a new buffer particle
                            !! nbuffer increases for total particle number increased
                            nbuffer = nbuffer + 1
                            k = getFixedPointIndex(point)
                            Particles(ntotal+nbuffer+nghost)%x = Particles(i)%x        &
                                - ( dot_product( (Particles(i)%x - fixedPoints(:, k)), &
                                                  normal(:, k) ) + dx * layer )        &
                                * normal(:, k)
                            Particles(ntotal+nbuffer+nghost)%State = 1
                            nghost = nghost + 1
                            Particles(ntotal+nbuffer+nghost)%x &
                                = calcGhostPosition(point_t(Particles(i)%x, 0))
                            Particles(ntotal+nbuffer+nghost)%State = 2
                            Particles(ntotal+nbuffer+nghost - 1)%neighborNum = 1
                            Particles(ntotal+nbuffer+nghost - 1)%neighborList(1) = ntotal+nbuffer+nghost
                        else !! Fluid Domain do not contain particle i
                            Particles(i)%State = -1
                            Particles(Particles(i)%neighborList(1))%State = -1
                        end if !! FluidDomain%contain(point)
                    end if !! .not. BufferDomainContain(point)
                end if !! Particles(i)%State == 1, Buffer Particle
            end do !! i
    
            ntotal = ntotal + nbuffer + nghost
        end if

    contains
        pure logical function BufferDomainContain(x) result(contain)
            type(point_t), intent(in) :: x

            contain = ( AllDomain%contain(x) ) .and. &
                      ( .not. (FluidDomain%contain(x)) )

        end function BufferDomainContain

        pure integer function getFixedPointIndex(x) result(index)
            type(point_t), intent(in) :: x
            real(8) :: distance, minDistance
            integer iter

            minDistance = huge(0._8)
            index = 0
            do iter = 1, 4
                distance = norm2(x%center - fixedPoints(:, i))
                if ( distance < minDistance ) then
                    minDistance = distance
                    index = iter
                end if
            end do

        end function getFixedPointIndex

        pure function calcGhostPosition(x) result(ghostPosition)
            type(point_t), intent(in) :: x
            real(8) :: ghostPosition(2)
            integer :: index

            index = getFixedPointIndex(x)
            ghostPosition = x%center                                     &
                + 2 * ( dot_product( (fixedPoints(:, index) - x%center), &
                                      normal(:, index)) )                &
                    * normal(:, index)

        end function calcGhostPosition

    end subroutine gen_undex_chamber_nrbc

    subroutine calc_undex_chamber_nrbc_property(Particles)
        type(Particle), intent(inout) :: Particles(:)
        integer :: ntotal
        type(point_t) :: point
        type(Particle) :: ghostParticle
        real(8), dimension(2, 4) :: fixedPoints, normal
        logical :: first_entry = .true.
        real(8) :: dx
        save fixedPoints, first_entry, dx

        integer i

        allocate(ghostParticle%x(2), ghostParticle%v(2))
        if ( first_entry ) then
            dx = abs(Particles(2)%x(1) - Particles(1)%x(1))

            !! Left
            fixedPoints(:, 1) = [real(8) :: -0.5, -0.5] &
                - [real(8) :: 0.5, 0]                   &
                - [real(8) :: dx / 2, 0]
            normal(:, 1) = [1, 0]

            !! Bottom
            fixedPoints(:, 2) = [real(8) :: -0.5, -0.5] &
                - [real(8) :: 0, 0.5]                   &
                - [real(8) :: 0, dx / 2]
            normal(:, 2) = [0, 1]

            !! Right
            fixedPoints(:, 3) = [real(8) :: -0.5, -0.5] &
                + [real(8) :: 0.5, 0]                   &
                + [real(8) :: dx / 2, 0]
            normal(:, 3) = [-1, 0]

            !! Top
            fixedPoints(:, 4) = [real(8) :: -0.5, -0.5] &
                + [real(8) :: 0, 0.5]                   &
                + [real(8) :: 0, dx / 2]
            normal(:, 4) = [0, -1]

            first_entry = .false.
        end if

        ntotal = size(Particles)
        do i = 1, ntotal
            if ( Particles(i)%State /= 1 ) cycle
            point = point_t(Particles(i)%x, 0)
        end do

    end subroutine calc_undex_chamber_nrbc_property

end module boundary_condition_m