module bc_m
    use ctrl_dict, only: Project, Config, Field
    use sph,       only: Particle, allocateParticleList
    use geometry_m
    implicit none

    private

    public :: gen_non_reflecting_bc, calc_nrbc_property
contains
    subroutine gen_non_reflecting_bc(ntotal, Particles)
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: Particles(:)

        select case (Project%nick)
        case ("undex_chamber")
            call gen_undex_chamber_nrbc(ntotal, Particles)
        end select

        call shrink(ntotal, Particles)

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
        type(Particle), allocatable :: Buffers(:), Ghosts(:)
        integer :: nbuffer, nghost
        type(rectangle_t) :: FluidDomain, AllDomain
        type(point_t) :: point
        real(8), dimension(2, 4) :: fixedPoints, normal
        integer :: entry_count = 0
        real(8) :: dx
        save FluidDomain, AllDomain, fixedPoints, entry_count, dx

        integer :: N, layer = 4
        integer i, k, l

        nbuffer = 0
        nghost = 0
        entry_count = entry_count + 1
        if ( entry_count == 1 ) then
            dx = abs(Particles(2)%x(1) - Particles(2)%x(2))
            FluidDomain = rectangle_t(   &
                [real(8) :: 0, 0], &
                [real(8) :: 1, 1], &
                0 )
            AllDomain = rectangle_t(     &
                [real(8) ::  0, 0], &
                ([real(8) :: 1, 1]  &
                    + 2 * layer * dx)  , &
                0 )

            !! Left
            fixedPoints(:, 1) = FluidDomain%center        &
                - [real(8) :: FluidDomain%length(1)/2, 0]
            normal(:, 1) = [1, 0]

            !! Bottom
            fixedPoints(:, 2) = FluidDomain%center        &
                - [real(8) :: 0, FluidDomain%length(2)/2]
            normal(:, 2) = [0, 1]

            !! Right
            fixedPoints(:, 3) = FluidDomain%center        &
                + [real(8) :: FluidDomain%length(1)/2, 0]
            normal(:, 3) = [-1, 0]

            !! Top
            fixedPoints(:, 4) = FluidDomain%center        &
                + [real(8) :: 0, FluidDomain%length(2)/2]
            normal(:, 4) = [0, -1]

            !!! This N is for particles' number in each layer
            N = int(1 / dx) + 1
            do l = 1, layer
                !! Left
                do i = 1, N
                    nbuffer = nbuffer + 1
                    Particles(ntotal+nbuffer)%x = [-0.5,  0.5] + [- l, 1 - i] * dx
                    Particles(ntotal+nbuffer)%State = 1
                    Particles(ntotal+nbuffer)%Type  = 6
                end do
                !! Bottom
                do i = 1, N
                    nbuffer = nbuffer + 1
                    Particles(ntotal+nbuffer)%x = [-0.5, -0.5] + [1 + i, - l] * dx
                    Particles(ntotal+nbuffer)%State = 1
                    Particles(ntotal+nbuffer)%Type  = 6
                end do
                !! Right
                do i = 1, N
                    nbuffer = nbuffer + 1
                    Particles(ntotal+nbuffer)%x = [ 0.5, -0.5] + [  l, 1 + i] * dx
                    Particles(ntotal+nbuffer)%State = 1
                    Particles(ntotal+nbuffer)%Type  = 6
                end do
                !! Top
                do i = 1, N
                    nbuffer = nbuffer + 1
                    Particles(ntotal+nbuffer)%x = [ 0.5,  0.5] + [1 - i,   l] * dx
                    Particles(ntotal+nbuffer)%State = 1
                    Particles(ntotal+nbuffer)%Type  = 6
                end do
            end do
            nghost = nbuffer
            call allocateParticleList(Ghosts, nghost, Field%dim, Field%pairNum)
            do i = ntotal + nbuffer + 1, ntotal + nbuffer + nghost
                Particles(i)%x = calcGhostPosition(point_t(Particles(i - nghost)%x, 0))
                Particles(i)%State = 2
                Particles(i - nghost)%neighborNum = 1
                Particles(i - nghost)%neighborList(1) = i
            end do
        else !! entry_count /= 1
            N = ntotal
            do i = 1, N
                point = point_t(Particles(i)%x, 0)
                if ( Particles(i)%State == 0 ) then !! Fluid Particle
                    if ( BufferDomainContain(point) ) then
                        !! Fluid Particle converts to Buffer Particle
                        Particles(i)%State = 1
                        nghost = nghost + 1
                        k = ntotal+nbuffer+nghost
                        Particles(k)%x = calcGhostPosition(point)
                        Particles(k)%State = 2
                        Particles(i)%neighborNum = 1
                        Particles(i)%neighborList(1) = k
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
                            k = ntotal+nbuffer+nghost
                            Particles(k)%x = newBufferParticlePosition(point)
                            Particles(k)%State = 1
                            Particles(k)%Type  = Particles(i)%Type
                            !! Add a new ghost particle
                            nghost = nghost + 1
                            k = ntotal+nbuffer+nghost
                            Particles(k)%x = calcGhostPosition(point)
                            Particles(k)%State = 2
                            Particles(k-1)%neighborNum = 1
                            Particles(k-1)%neighborList = 0
                            Particles(k-1)%neighborList(1) = k
                        else !! Fluid Domain do not contain particle i
                            Particles(i)%State = -1
                            Particles(Particles(i)%neighborList(1))%State = -1
                        end if !! FluidDomain%contain(point)
                    end if !! .not. BufferDomainContain(point)
                end if !! Particles(i)%State == 1, Buffer Particle
            end do !! i

        end if !! entry_count == 1

        ntotal = ntotal + nbuffer + nghost

    contains
        pure logical function BufferDomainContain(P) result(contain)
            type(point_t), intent(in) :: P

            contain = ( AllDomain%contain(P) ) .and. &
                      ( .not. (FluidDomain%contain(P)) )

        end function BufferDomainContain

        pure integer function getFixedPointIndex(P) result(index)
            type(point_t), intent(in) :: P
            real(8) :: distance, minDistance
            integer iter

            minDistance = huge(0._8)
            index = 0
            do iter = 1, 4
                distance = norm2(P%center - fixedPoints(:, iter))
                if ( distance < minDistance ) then
                    minDistance = distance
                    index = iter
                end if
            end do

        end function getFixedPointIndex

        pure function newBufferParticlePosition(P) result(bufferPosition)
            type(point_t), intent(in) :: P
            real(8) :: bufferPosition(Field%dim)
            integer :: index

            index = getFixedPointIndex(P)
            bufferPosition = P%center                                &
                - ( dot_product( (P%center - fixedPoints(:, index)), &
                                  normal(:, index) ) + dx * layer )  &
                * normal(:, index)

        end function newBufferParticlePosition

        pure function calcGhostPosition(P) result(ghostPosition)
            type(point_t), intent(in) :: P
            real(8) :: ghostPosition(Field%dim)
            integer :: index

            index = getFixedPointIndex(P)
            ghostPosition = P%center                                     &
                + 2 * ( dot_product( (fixedPoints(:, index) - P%center), &
                                      normal(:, index)) )                &
                    * normal(:, index)

        end function calcGhostPosition

    end subroutine gen_undex_chamber_nrbc

    subroutine calc_undex_chamber_nrbc_property(Particles)
        use tools_m, only: dyadic_product
        type(Particle), intent(inout) :: Particles(:)
        integer :: ntotal
        real(8), allocatable :: v(:, :)

        integer i, d, g

        ntotal = size(Particles)
        allocate(v(Field%dim, ntotal))
        do i = 1, ntotal
            v(:, i) = Particles(i)%v(:)
        end do
        do i = 1, ntotal
            if ( Particles(i)%State /= 1 ) cycle
            g = Particles(i)%neighborList(1)
            do d = 1, Field%dim
                Particles(i)%v(d) = solveBufferProperty(Particles(i)%x, Particles(g), v(d, :))
            end do
            Particles(i)%Mass = Particles(1)%Mass
            Particles(i)%Density = solveBufferProperty(Particles(i)%x, Particles(g), Particles(:)%Density)
            Particles(i)%Pressure = solveBufferProperty(Particles(i)%x, Particles(g), Particles(:)%Pressure)
            Particles(i)%InternalEnergy = solveBufferProperty(Particles(i)%x, Particles(g), Particles(:)%InternalEnergy)
            Particles(i)%SoundSpeed = solveBufferProperty(Particles(i)%x, Particles(g), Particles(:)%SoundSpeed)
            Particles(i)%divergenceVelocity = solveBufferProperty(Particles(i)%x, Particles(g), Particles(:)%divergenceVelocity)
        end do

    contains
        function solveBufferProperty(coor, ghost, Property) result(bufferProperty)
            real(8), intent(in) :: coor(:)
            type(Particle), intent(in) :: ghost
            real(8), intent(in) :: Property(:)
            real(8), allocatable :: A(:,:), Solve(:), h(:)
            real(8) :: bufferProperty
            integer :: m, flag

            integer j, l

            m = Field%dim + 1
            allocate(A(m, m), Solve(m), h(m), source=0._8)
            do j = 1, ghost%neighborNum
                l = ghost%neighborList(j)
                associate (aux => [ghost%w(j), ghost%dwdx(:, j)] &
                    * Particles(l)%Mass/Particles(l)%Density)
                A = A + dyadic_product(aux, [1._8, Particles(l)%x-ghost%x])
                Solve = Solve + aux * Property(l)
                end associate
            end do

            call DGESV(m,1, A,m, h, Solve,m, flag)

            bufferProperty = Solve(1) + dot_product((coor - ghost%x), Solve(2:m))
            
        end function solveBufferProperty

    end subroutine calc_undex_chamber_nrbc_property

    subroutine shrink(ntotal, Particles)
        integer, intent(inout) :: ntotal
        type(Particle), intent(inout) :: Particles(:)
        type(Particle) :: P

        integer i, j

        i = 1
        j = ntotal
        do while( i <= j )
            do while( i <= j .and. Particles(i)%State /= -1 )
                i = i + 1
            end do
            do while( i <= j .and. Particles(j)%State == -1 )
                j = j - 1
            end do
            if ( i < j ) then
                P = Particles(i)
                Particles(i) = Particles(j)
                Particles(j) = P
                i = i + 1
                j = j - 1
            end if
        end do

        ntotal = j
        
    end subroutine shrink

end module bc_m