module bc_m
    use ctrl_dict, only: Project, Config, Field
    use sph,       only: Particle, allocateParticleList
    use geometry_m
    implicit none
    private

    type Buffer_t
        type(Particle) :: b !! Buffer
        type(Particle) :: g !! Ghost
    end type Buffer_t

    public :: gen_non_reflecting_bc
contains
    subroutine gen_non_reflecting_bc(ntotal, Particles, nbuffer, Buffers)
        integer, intent(inout) :: ntotal, nbuffer
        type(Particle), intent(inout) :: Particles(:), Buffers(:)

        select case (Project%nick)
        case ("undex_chamber")
            call gen_undex_chamber_nrbc(ntotal, Particles, nbuffer, Buffers)
        end select

    end subroutine gen_non_reflecting_bc

    !!! Open Boundary Non-Reflecting Condition
    subroutine gen_undex_chamber_nrbc(ntotal, Particles, nbuffer, Buffers)
        use APS_M, only: BGGS
        integer, intent(inout) :: ntotal, nbuffer
        type(Particle), intent(inout) :: Particles(:), Buffers(:)
        type(Particle), allocatable :: Ghosts(:)
        integer :: N
        type(rectangle_t) :: FluidDomain, AllDomain
        real(8), dimension(2, 4) :: fixedPoints, normal
        integer :: entry_count = 0
        real(8) :: dx
        save FluidDomain, AllDomain, fixedPoints, entry_count, dx

        integer :: layer = 4
        integer i, l

        entry_count = entry_count + 1

        if ( entry_count == 1 ) then
            nbuffer = 0
            dx = abs(Particles(2)%x(1) - Particles(2)%x(2))
            FluidDomain = rectangle_t( &
                [real(8) :: 0, 0],     &
                [real(8) :: 1, 1],     &
                0 )
            AllDomain = rectangle_t(   &
                [real(8) ::  0, 0],    &
                ([real(8) :: 1, 1]     &
                    + 2 * layer * dx), &
                0 )

            !! Left
            fixedPoints(:, 1) = FluidDomain%center &
                - [real(8) :: FluidDomain%length(1)/2, 0]
            normal(:, 1) = [1, 0]

            !! Bottom
            fixedPoints(:, 2) = FluidDomain%center &
                - [real(8) :: 0, FluidDomain%length(2)/2]
            normal(:, 2) = [0, 1]

            !! Right
            fixedPoints(:, 3) = FluidDomain%center &
                + [real(8) :: FluidDomain%length(1)/2, 0]
            normal(:, 3) = [-1, 0]

            !! Top
            fixedPoints(:, 4) = FluidDomain%center &
                + [real(8) :: 0, FluidDomain%length(2)/2]
            normal(:, 4) = [0, -1]

            FluidDomain%length = FluidDomain%length + 0.1 * dx

            !!! This N is for particles' number in each layer
            N = int(1 / dx) + 1
            do l = 1, layer
                !! Left
                do i = 1, N
                    nbuffer = nbuffer + 1
                    Buffers(nbuffer)%x = [-0.5,  0.5] + [- l, 1 - i] * dx
                    Buffers(nbuffer)%State = 1
                    Buffers(nbuffer)%Type  = 6
                end do
                !! Bottom
                do i = 1, N
                    nbuffer = nbuffer + 1
                    Buffers(nbuffer)%x = [-0.5, -0.5] + [i - 1, - l] * dx
                    Buffers(nbuffer)%State = 1
                    Buffers(nbuffer)%Type  = 6
                end do
                !! Right
                do i = 1, N
                    nbuffer = nbuffer + 1
                    Buffers(nbuffer)%x = [ 0.5, -0.5] + [  l, i - 1] * dx
                    Buffers(nbuffer)%State = 1
                    Buffers(nbuffer)%Type  = 6
                end do
                !! Top
                do i = 1, N
                    nbuffer = nbuffer + 1
                    Buffers(nbuffer)%x = [ 0.5,  0.5] + [1 - i,   l] * dx
                    Buffers(nbuffer)%State = 1
                    Buffers(nbuffer)%Type  = 6
                end do
            end do !! l
            do i = 1, nbuffer
                Ghosts(i)%x = calcGhostPosition(Buffers(i))
            end do
        else !! entry_count /= 1
        ! end if !! entry_count == 1

            do i = 1, ntotal !! Fluid Particle
                if ( BufferDomainContain(Particles(i)) ) then
                    !! Fluid Particle converts to Buffer Particle
                    Particles(i)%State = -1
                    nbuffer = nbuffer + 1
                    Buffers(nbuffer)%x = Particles(i)%x
                    Buffers(nbuffer)%State = 0
                    Buffers(nbuffer)%Type = Particles(i)%Type
                    Ghosts(nbuffer)%x = calcGhostPosition(Buffers(nbuffer))
                end if
            end do

            N = nbuffer
            do i = 1, N !! Buffer Particle
                if ( .not. BufferDomainContain(Buffers(i)) ) then
                    if ( FluidDomain%contain(point_t(Buffers(i)%x, 0)) ) then
                        !! Buffer Particle converts to Fluid Particle
                        Buffers(i)%State = -1
                        ntotal = ntotal + 1
                        Particles(ntotal) = Buffers(i)

                        !! Add a new buffer particle
                        !! nbuffer increases for total particle number increased
                        nbuffer = nbuffer + 1
                        Buffers(nbuffer)%x = newBufferParticlePosition(Buffers(i))
                        Buffers(nbuffer)%State = 0
                        Buffers(nbuffer)%Type  = Particles(ntotal)%Type
                        !! Add a new ghost particle
                        Ghosts(nbuffer)%x = calcGhostPosition(Buffers(nbuffer))
                    else !! Fluid Domain do not contain particle i
                        Buffers(i)%State = -1
                    end if !! FluidDomain%contain(point)
                end if
            end do

        end if !! entry_count == 1

        ! call shrink(Particles, ntotal)
        ! call bufferShrink(Buffers, nbuffer)
        ! Buffers = [Buffers(1:nbuffer)]

        call allocateParticleList(Ghosts, nbuffer, Field%dim, Field%pairNum)
        do i = 1, nbuffer
            Ghosts(i)%SmoothingLength = 1.5 * dx
        end do
        call BGGS(Ghosts, Particles(1:ntotal))
        do i = 1, nbuffer
            Ghosts(i) = Ghosts(i)
        end do

        do i = 1, nbuffer
            call solveBufferProperty(Buffers(i), Ghosts(i), Particles)
            Particles(ntotal + i) = Buffers(i)
        end do

        N = nbuffer

    contains
        pure logical function BufferDomainContain(P) result(contain)
            type(Particle), intent(in) :: P
            type(point_t) :: point
            integer :: err

            allocate(point%center(Field%dim), source=0._8, stat=err)
            if ( err /= 0 ) error stop "BufferDomainContain: allocate point%center failed"
            point = point_t(P%x, 0)
            contain = ( AllDomain%contain(point) ) .and. &
                      ( .not. (FluidDomain%contain(point)) )
            deallocate(point%center)

        end function BufferDomainContain

        pure integer function getFixedPointIndex(P) result(index)
            type(Particle), intent(in) :: P
            real(8) :: distance, minDistance
            integer iter

            minDistance = huge(0._8)
            index = 0
            do iter = 1, 4
                distance = norm2(P%x - fixedPoints(:, iter))
                if ( distance < minDistance ) then
                    minDistance = distance
                    index = iter
                end if
            end do

        end function getFixedPointIndex

        pure function newBufferParticlePosition(P) result(bufferPosition)
            type(Particle), intent(in) :: P
            real(8) :: bufferPosition(Field%dim)
            integer :: index

            index = getFixedPointIndex(P)
            bufferPosition = P%x                                    &
                - ( dot_product( (P%x - fixedPoints(:, index)),     &
                                  normal(:, index) ) + dx * layer ) &
                * normal(:, index)

        end function newBufferParticlePosition

        pure function calcGhostPosition(P) result(ghostPosition)
            type(Particle), intent(in) :: P
            real(8) :: ghostPosition(Field%dim)
            integer :: index

            index = getFixedPointIndex(P)
            ghostPosition = P%x                                     &
                + 2 * ( dot_product( (fixedPoints(:, index) - P%x), &
                                      normal(:, index)) )           &
                    * normal(:, index)

        end function calcGhostPosition

    end subroutine gen_undex_chamber_nrbc

    subroutine solveBufferProperty(buffer, ghost, Particles)
        use tools_m,  only: dyadic_product
        ! use kernel_m, only: kernel
        type(Particle), intent(inout) :: buffer
        type(Particle), intent(in)    :: ghost
        type(Particle), intent(in)    :: Particles(:)
        real(8), allocatable :: A(:,:), Solve(:, :), h(:)
        ! real(8) :: wi
        integer :: m, flag

        integer j, l, d

        m = Field%dim + 1
        allocate(A(m, m), Solve(m, Field%dim+6), h(m), source=0._8)
        h = 0
        do j = 1, ghost%neighborNum
            l = ghost%neighborList(j)
            associate (aux => [ghost%w(j), ghost%dwdx(:, j)] &
                * Particles(l)%Mass/Particles(l)%Density)
            h = h + aux * [1._8, Particles(l)%x-ghost%x]
            A = A + dyadic_product(aux, [1._8, Particles(l)%x-ghost%x])
            Solve(:, 1) = Solve(:, 1) + aux * Particles(l)%Mass
            Solve(:, 2) = Solve(:, 2) + aux * Particles(l)%Density
            Solve(:, 3) = Solve(:, 3) + aux * Particles(l)%Pressure
            Solve(:, 4) = Solve(:, 4) + aux * Particles(l)%InternalEnergy
            Solve(:, 5) = Solve(:, 5) + aux * Particles(l)%SoundSpeed
            Solve(:, 6) = Solve(:, 6) + aux * Particles(l)%divergenceVelocity
            do d = 7, Field%dim + 6
                Solve(:, d) = Solve(:, d) + aux * Particles(l)%v(d)
            end do
            end associate
        end do

        Solve(1, :) = Solve(1, :) / h(1)
        do d = 2, Field%dim + 1
            Solve(d, :) = Solve(d, :) / h(d)
        end do

        h = 0
        do j = 1, Field%dim + 6
            call DGESV(m,1, A,m, h, Solve(:, j),m, flag)
        end do

        buffer%Mass               = Solve(1, 1) + dot_product((buffer%x - ghost%x), Solve(2:m, 1))
        buffer%Density            = Solve(1, 2) + dot_product((buffer%x - ghost%x), Solve(2:m, 2))
        buffer%Pressure           = Solve(1, 3) + dot_product((buffer%x - ghost%x), Solve(2:m, 3))
        buffer%InternalEnergy     = Solve(1, 4) + dot_product((buffer%x - ghost%x), Solve(2:m, 4))
        buffer%SoundSpeed         = Solve(1, 5) + dot_product((buffer%x - ghost%x), Solve(2:m, 5))
        buffer%divergenceVelocity = Solve(1, 6) + dot_product((buffer%x - ghost%x), Solve(2:m, 6))
        do d = 1, Field%dim
            buffer%v(d) = Solve(1, d + 6) + dot_product((buffer%x - ghost%x), Solve(2:m, d + 6))
        end do

    end subroutine solveBufferProperty

    subroutine shrink(Particles, ntotal)
        type(Particle), intent(inout) :: Particles(:)
        integer, intent(inout) :: ntotal
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

    ! subroutine bufferShrink(Particles, ntotal)
    !     type(Buffer_t), intent(inout) :: Particles(:)
    !     integer, intent(inout) :: ntotal
    !     type(Buffer_t) :: P

    !     integer i, j

    !     i = 1
    !     j = ntotal
    !     do while ( i <= j )
    !         do while( i <= j .and. Particles(i)%State /= -1 )
    !             i = i + 1
    !         end do
    !         do while( i <= j .and. Particles(j)%State == -1 )
    !             j = j - 1
    !         end do
    !         if ( i < j ) then
    !             P = Particles(i)
    !             Particles(i) = Particles(j)
    !             Particles(j) = P
    !             i = i + 1
    !             j = j - 1
    !         end if
    !     end do

    !     ntotal = j

    ! end subroutine bufferShrink

end module bc_m