module BGGS_m
    use ctrl_dict, only: Field, Config
    use SPH,       only: Particle
    implicit none
    private

    type cell_t
        integer :: numParticles
        integer, allocatable :: particleList(:)
    end type cell_t

    type grid_t
        type(cell_t), allocatable :: grid(:, :, :)
        integer, allocatable :: cellNum(:)             !! Number of grid cells
        real(8), allocatable :: maxCoor(:), minCoor(:) !! Maximum and minium grid cell coordinates
    contains
        procedure :: initialize => initialize
        procedure :: locate     => locate
    end type grid_t

    real(8), parameter :: GMR = 1.2 !! Grid Magnification Ratio
    integer, parameter :: NPG = 3   !! Number of particles Per Grid cell

    public :: BGGS

contains
    pure function locate(self, p) result(cell)
        use tools_m, only: to_string
        class(grid_t),  intent(in) :: self
        type(Particle), intent(in) :: p
        integer, allocatable :: cell(:)

        integer d

        allocate( cell(Field%dim), source = 0)

        do d = 1, Field%dim
            if ((p%x(d) > self%maxCoor(d)) .or. (p%x(d) < self%minCoor(d))) then
                error stop "particle out of range"
            else
                cell(d) = int( ( (p%x(d) - self%minCoor(d))       &
                            / (self%maxCoor(d)-self%minCoor(d)) ) &
                            * self%cellNum(d) + 1 )
            end if
        end do

    end function locate

    subroutine initialize(self, Particles)
        class(grid_t), intent(inout) :: self
        type(Particle), intent(in) :: Particles(:)
        integer, allocatable :: cell(:)
        integer :: ntotal
        logical :: firstEntry = .true.

        integer i, d, dd, ddd, n, maxn
        save firstEntry, n

        ntotal = size(Particles)
        allocate( self%maxCoor(Field%dim), source =-huge(0._8) )
        allocate( self%minCoor(Field%dim), source = huge(0._8) )
        allocate( self%cellNum(Field%dim), source = 0 )
        allocate( cell(Field%dim), source = 0 )

        do i = 1, ntotal
            do d = 1, Field%dim
                if ( Particles(i)%x(d) > self%maxCoor(d) ) self%maxCoor(d) = Particles(i)%x(d)
                if ( Particles(i)%x(d) < self%minCoor(d) ) self%minCoor(d) = Particles(i)%x(d)
            end do
        end do

        associate( delta => (self%maxCoor-self%minCoor) )
            self%maxCoor = self%maxCoor + (GMR - 1) * delta
            self%minCoor = self%minCoor - (GMR - 1) * delta
            !!! number of grid cells in x-, y- and z-direction:
            select case (Field%dim)
            case (1)
                self%cellNum(1) = min((ntotal)/NPG + 1, 1000)
                allocate( self%grid(self%cellNum(1), 1, 1 ) )
                if ( firstEntry ) n = (ntotal / product(self%cellNum)) + 1
                do d = 1, self%cellNum(1)
                    self%grid(d, 1, 1)%numParticles = 0
                    allocate( self%grid(d, 1, 1)%particleList(n), source = 0 )
                end do
                do i = 1, ntotal
                    cell = self%locate(Particles(i))
                    self%grid(cell(1), 1, 1)%numParticles &
                        = self%grid(cell(1), 1, 1)%numParticles + 1
                    if ( self%grid(cell(1), 1, 1)%numParticles <= n ) then
                        self%grid(cell(1), 1, 1)% &
                            particleList(self%grid(cell(1), 1, 1)%numParticles) = i
                    else
                        self%grid(cell(1), 1, 1)%particleList &
                            = [self%grid(cell(1), 1, 1)%particleList, i]
                    end if
                    if ( self%grid(cell(1), 1, 1)%numParticles > n ) then
                        maxn = self%grid(cell(1), 1, 1)%numParticles
                    end if
                end do
            case (2)
                self%cellNum(1) = min(int(((ntotal) * delta(1) &
                                / (delta(2)*NPG))**(1._8/2) ) + 1, 1000)
                self%cellNum(2) = min(int(self%cellNum(1)*delta(2)/delta(1)) + 1, 1000)
                allocate( self%grid(self%cellNum(1), self%cellNum(2), 1 ) )
                if ( firstEntry ) n = (ntotal / product(self%cellNum)) + 1
                do d = 1, self%cellNum(1)
                    do dd = 1, self%cellNum(2)
                        self%grid(d, dd, 1)%numParticles = 0
                        allocate( self%grid(d, dd, 1)%particleList(n), source = 0 )
                    end do
                end do
                do i = 1, ntotal
                    cell = self%locate(Particles(i))
                    self%grid(cell(1), cell(2), 1)%numParticles &
                        = self%grid(cell(1), cell(2), 1)%numParticles + 1
                    if ( self%grid(cell(1), cell(2), 1)%numParticles <= n ) then
                        self%grid(cell(1), cell(2), 1)% &
                            particleList(self%grid(cell(1), cell(2), 1)%numParticles) = i
                    else
                        self%grid(cell(1), cell(2), 1)%particleList &
                            = [self%grid(cell(1), cell(2), 1)%particleList, i]
                    end if
                    if ( self%grid(cell(1), cell(2), 1)%numParticles > n ) then
                        maxn = self%grid(cell(1), cell(2), 1)%numParticles
                    end if
                end do
            case (3)
                self%cellNum(1) = min(int(((ntotal) * delta(1) * delta(1) &
                                / (delta(2) * delta(3) * NPG))**(1._8/3))+ 1, 1000)
                self%cellNum(2) = min(int(self%cellNum(1) * delta(2) / delta(1)) + 1, 1000)
                self%cellNum(3) = min(int(self%cellNum(1) * delta(3) / delta(1)) + 1, 1000)
                allocate( self%grid(self%cellNum(1), self%cellNum(2), self%cellNum(3)) )
                if ( firstEntry ) n = (ntotal / product(self%cellNum)) + 1
                do d = 1, self%cellNum(1)
                    do dd = 1, self%cellNum(2)
                        do ddd = 1, self%cellNum(3)
                            self%grid(d, dd, ddd)%numParticles = 0
                            allocate( self%grid(d, dd, ddd)%particleList(n), source = 0 )
                        end do
                    end do
                end do
                do i = 1, ntotal
                    cell = self%locate(Particles(i))
                    self%grid(cell(1), cell(2), cell(3))%numParticles &
                        = self%grid(cell(1), cell(2), cell(3))%numParticles + 1
                    if ( self%grid(cell(1), cell(2), cell(3))%numParticles <= n ) then
                        self%grid(cell(1), cell(2), cell(3))% &
                            particleList(self%grid(cell(1), cell(2), cell(3))%numParticles) = i
                    else
                        self%grid(cell(1), cell(2), cell(3))%particleList &
                            = [self%grid(cell(1), cell(2), cell(3))%particleList, i]
                    end if
                    if ( self%grid(cell(1), cell(2), cell(3))%numParticles > n ) then
                        maxn = self%grid(cell(1), cell(2), cell(3))%numParticles
                    end if
                end do
            end select
        end associate

        n = maxn

    end subroutine initialize

    subroutine BGGS(Targets, Sources, skipItsSelf)
        use kernel_m, only: kernel
        type(Particle), intent(inout) :: Targets(:)
        type(Particle), intent(in)    :: Sources(:)
        logical, intent(in), optional :: skipItsSelf
        type(grid_t) :: Grid
        integer :: numTargets, numPairs
        integer :: scale_k
        integer :: maxCell(3), minCell(3)
        integer, allocatable :: cell(:), cellNumPerHSML(:)
        real(8), allocatable :: dx(:)
        real(8) :: dr, r, mhsml

        integer i, j, k, d, xcell, ycell, zcell

        select case (Config%skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        case default
            scale_k = 1
        end select

        numTargets = size(Targets)
        numPairs   = size(Targets(1)%neighborList)

        allocate( cell(Field%dim), source = 0 )
        allocate( cellNumPerHSML(Field%dim), source = 0 )
        allocate( dx(Field%dim), source = 0._8 )

        call Grid%initialize(Sources)

        !$OMP PARALLEL DO PRIVATE(i, j, k, d, xcell, ycell, zcell) &
        !$OMP PRIVATE(cell, maxCell, minCell, cellNumPerHSML)      &
        !$OMP PRIVATE(dx, dr, r, mhsml)                            &
        !$OMP SCHEDULE(dynamic, Config%chunkSize)
        do i = 1, numTargets
            cell = Grid%locate(Targets(i))
            cellNumPerHSML = int( Targets(i)%SmoothingLength &
                / ((Grid%maxCoor - Grid%minCoor)             &
                    / Grid%cellNum) + 1 )
            maxCell = 1
            minCell = 1
            do d = 1, Field%dim
                maxCell(d) = min(cell(d) + cellNumPerHSML(d), Grid%cellNum(d))
                minCell(d) = max(cell(d) - cellNumPerHSML(d), 1)
            end do
            do zcell = minCell(3), maxCell(3)
                do ycell = minCell(2), maxCell(2)
                    do xcell = minCell(1), maxCell(1)
                        do k = 1, Grid%grid(xcell, ycell, zcell)%numParticles
                            j = Grid%grid(xcell, ycell, zcell)%particleList(k)
                            if ( present(skipItsSelf) ) then
                                if ( skipItsSelf .and. (i == j) ) cycle
                            end if ! present(skipItsSelf)
                            dx = Targets(i)%x(1) - Sources(j)%x(1)
                            dr = dx(1) ** 2
                            do d = 2, Field%dim
                                dx(d) = Targets(i)%x(d) - Sources(j)%x(d)
                                dr = dr + dx(d) ** 2
                            end do !! d
                            r = dr ** 0.5
                            mhsml = 0.5 * (Targets(i)%SmoothingLength &
                                         + Sources(j)%SmoothingLength)
                            if ( r < mhsml * scale_k ) then
                                Targets(i)%neighborNum = Targets(i)%neighborNum + 1
                                if ( Targets(i)%neighborNum <= numPairs ) then
                                    Targets(i)%neighborList(Targets(i)%neighborNum) = j
                                else !! Targets(i)%neighborNum > numPairs
                                    Targets(i)%neighborList = [Targets(i)%neighborList, j]
                                end if !! Targets(i)%neighborNum <= numPairs
                                call kernel(r, dx, mhsml,                 &
                                    Targets(i)%w(Targets(i)%neighborNum), &
                                    Targets(i)%dwdx(:, Targets(i)%neighborNum))
                            end if !! r < mhsml * scale_k
                        end do !! k
                    end do !! xcell
                end do !! ycell
            end do !! zcell
        end do !! i

    end subroutine BGGS

end module BGGS_m