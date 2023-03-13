#include "../macro.h"
module nnps_m
    use, intrinsic :: iso_fortran_env, only: err => error_unit
    use ctrl_dict, only: dim, skf, chunkSize
    use sph
    use kernel_m,  only: kernel
    use tools_m,   only: to_string, round, print_error, print_warning, ESC
    implicit none

    private

    public :: search_particles, print_statistics

contains
    subroutine search_particles(nnps, Particles)
        integer, intent(in)  :: nnps
        type(Particle), intent(inout) :: Particles(:)

#if CHECK_NEIGHBOR_LIST
        integer :: i, j, k
#endif

        select case (nnps)
        case (1)
            call direct_search(Particles)
        case (2)
            call link_list_search(Particles)
        case (3)
            call tree_search(Particles)
        end select

#if CHECK_NEIGHBOR_LIST
        do i = 1, size(x, 2)
            do k = 1, Particles(i)%neighborNum
                j = Particles(i)%neighborList(k)
                if ( j == 0 ) then
                    write(err,"(A)",advance="no") ESC//"[31m"
                    write(err,"(2(A,I0),A)") "Neighbor ", k, &
                                            " of Particle ", i, " is 0."
                    write(err,"(A)",advance="no") ESC//"[0m"
                end if
            end do
        end do
#endif

    end subroutine search_particles

    subroutine direct_search(P)
        !$ use omp_lib
        ! include "mpif.h"
        type(Particle), intent(inout) :: P(:)
        integer :: ntotal, kpair
        real(8) :: this_w, this_dwdx(dim)
        integer :: scale_k
        real(8) :: dx(dim), dr, r, mhsml

        integer i, j, d !! Loop variables

        ! integer :: tid !! OMP Parallel parameters
        ! integer :: nprocs, pid, ierror, remain  !! MPI Parallel parameters

        select case (skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        case default
            scale_k = 1
        end select

        ntotal = size(P)
        kpair = size(P(1)%neighborList)

        do i = 1, ntotal - 1
            do j = i + 1, ntotal
                ! dx = x(:, i) - x(:, j)
                ! r  = norm2(dx)
                dx(1) = P(i)%x(1) - P(j)%x(1)
                dr = dx(1)*dx(1)
                do d = 2, dim
                    dx(d) = P(i)%x(d) - P(j)%x(d)
                    dr = dr + dx(d)*dx(d)
                end do
                r = sqrt(dr)
                mhsml = 0.5_8 * ( P(i)%SmoothingLength + P(j)%SmoothingLength )
                if ( r < mhsml * scale_k ) then
                    !!! Neighbor particle ant tottal neighbor number for each particle
                    P(i)%neighborNum = P(i)%neighborNum + 1
                    P(j)%neighborNum = P(j)%neighborNum + 1
#if CHECK_NEIGHBOR_NUM
                    if ( P(i)%neighborNum > kpair .or. P(j)%neighborNum > kpair ) then
                        error stop "Too many neighbors"
                    end if
#endif
                    P(i)%neighborList(P(i)%neighborNum) = j
                    P(j)%neighborList(P(j)%neighborNum) = i
                    !!! Kernel and derivations of kernel
                    call kernel(r, dx, mhsml, this_w, this_dwdx)
                    P(i)%w(P(i)%neighborNum) = this_w
                    P(j)%w(P(j)%neighborNum) = this_w
                    P(i)%dwdx(:, P(i)%neighborNum) = this_dwdx
                    P(j)%dwdx(:, P(j)%neighborNum) = -this_dwdx
                end if
            end do
        end do

    end subroutine direct_search

    subroutine link_list_search(P)
        type(Particle), intent(inout) :: P(:)
        integer :: ntotal, kpair
        real(8) :: this_w, this_dwdx(dim)
        real(8) :: maxhsml, mhsml
        integer, allocatable   :: grid(:, :, :)
        integer, allocatable   :: cell_index(:, :), cell_data(:)
        integer :: cell(3), xcell, ycell, zcell, &
                   minxcell(3), maxxcell(3), &
                   cell_num(dim), ghsmlx(dim)
        real(8) :: dr, r, dx(dim), grid_max_coor(dim), grid_min_coor(dim)
        real(8) :: scale
        integer :: i, j, d, scale_k

        select case (skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        case default
            scale_k = 1
        end select

        ntotal = size(P)
        kpair = size(P(1)%neighborList)

        allocate(cell_index(3, ntotal), source=0)
        allocate(cell_data(ntotal), source=0)

        grid_max_coor(:) = -huge(0._8)
        grid_min_coor(:) =  huge(0._8)
        do i = 1, ntotal
            do d = 1, dim
                if ( P(i)%x(d) > grid_max_coor(d) ) grid_max_coor(d) = P(i)%x(d)
                if ( P(i)%x(d) < grid_min_coor(d) ) grid_min_coor(d) = P(i)%x(d)
            end do
        end do
        scale = 1.2
        associate( d => (grid_max_coor-grid_min_coor) )
            grid_max_coor = grid_max_coor + (scale - 1) * d
            grid_min_coor = grid_min_coor - (scale - 1) * d
        end associate
        
        !!! initialize grid:
        call initialize_grid(ntotal, maxval(P%SmoothingLength), grid, cell_num, ghsmlx, grid_max_coor, grid_min_coor)

        !!! position particles on grid and create linked list:
        do i = 1, ntotal
            call grid_geometry(i, P(i)%x(:), cell_num, grid_max_coor, grid_min_coor, cell)
            cell_index(:, i) = cell(:)
            cell_data(i) = grid(cell(1), cell(2), cell(3))
            grid(cell(1), cell(2), cell(3)) = i
        end do

        !$OMP PARALLEL DO PRIVATE(i, j, d, cell, xcell, ycell, zcell, minxcell, maxxcell, dx, dr, r, mhsml, this_w, this_dwdx) &
        !$OMP SHARED(grid, cell_index, cell_data, P, ghsmlx, cell_num) &
        !$OMP SCHEDULE(DYNAMIC, chunkSize)
        !!! determine interaction parameters:
        do i = 1, ntotal - 1 !! loop over all particles but the last one
            !!! determine range of grid to go through:
            maxxcell(:) = 1
            minxcell(:) = 1
            do d = 1, dim
                maxxcell(d) = min(cell_index(d, i) + ghsmlx(d), cell_num(d))
                minxcell(d) = max(cell_index(d, i) - ghsmlx(d), 1)
            end do

            !!! search grid:
            do zcell = minxcell(3), maxxcell(3) !! Loop over all cells at z direction
                do ycell = minxcell(2), maxxcell(2) !! Loop over all cells at y direction
                    do xcell = minxcell(1), maxxcell(1) !! Loop over all cells at x direction
                        j = grid(xcell, ycell, zcell) !! Fetch Particle j from grid
                        do while (j > i)
                            !!! Calculate distance between particle i and j
                            dx(1) = P(i)%x(1) - P(j)%x(1)
                            dr = dx(1)*dx(1)
                            do d = 2, dim
                                dx(d) = P(i)%x(d) - P(j)%x(d)
                                dr = dr + dx(d)*dx(d)
                            end do
                            r = sqrt(dr)
                            mhsml = 0.5_8 * (P(i)%SmoothingLength + P(j)%SmoothingLength)
                            if (r < mhsml * scale_k) then
                                !!! Neighbor particle ant tottal neighbor number for each particle
                                !$OMP CRITICAL
                                P(i)%neighborNum = P(i)%neighborNum + 1
                                P(j)%neighborNum = P(j)%neighborNum + 1
#if CHECK_NEIGHBOR_NUM
                                if ( P(i)%neighborNum > kpair .or. P(j)%neighborNum > kpair ) then
                                    write(*,*) i, j, P([i,j])%neighborNum, kpair
                                    error stop "Too many neighbors"
                                end if
#endif
                                P(i)%neighborList(P(i)%neighborNum) = j
                                P(j)%neighborList(P(j)%neighborNum) = i
                                !!! Kernel and derivations of kernel
                                call kernel(r, dx, mhsml, this_w, this_dwdx)
                                P(i)%w(P(i)%neighborNum) = this_w
                                P(j)%w(P(j)%neighborNum) = this_w
                                P(i)%dwdx(:, P(i)%neighborNum) = this_dwdx
                                P(j)%dwdx(:, P(j)%neighborNum) = -this_dwdx
                                !$OMP END CRITICAL
                            end if !! r < mhsml * scale_k   
                            j = cell_data(j)
                        end do !! j
                    end do !! xcell
                end do !! ycell
            end do !! zcell
        end do !! i
        !$OMP END PARALLEL DO

        deallocate(grid)
        deallocate(cell_index, cell_data)

    end subroutine link_list_search

    pure subroutine initialize_grid(ntotal, hsml, grid, cell_num, ghsmlx, grid_max_coor, grid_min_coor)
        integer, intent(in) :: ntotal                                !! Number of particles
        real(8), intent(in) :: hsml                                  !! Smoothing length
        integer, intent(inout), allocatable :: grid(:, :, :)         !! Grid
        integer, intent(inout) :: cell_num(:)                        !! Number of grid cells
        integer, intent(inout) :: ghsmlx(:)                          !! Smoothing length
        real(8), intent(inout) :: grid_max_coor(:), grid_min_coor(:) !! Maximum and minium grid cell coordinates
        integer, parameter     :: nppg = 3  !! averaged number of particles per grid cell

        !!! initialize parameters: maximum number of grid cells

        !!! range of sorting grid

        associate( d=>(grid_max_coor(:) - grid_min_coor(:)) )
            !!! number of grid cells in x-, y- and z-direction:
            if (dim == 1) then
                cell_num(1) = min((ntotal)/nppg + 1, 1000)
            else if (dim == 2) then
                cell_num(1) = min(int(((ntotal)*d(1) &
                                /(d(2)*nppg))**(1._8/2) )     + 1, 1000)
                cell_num(2) = min(int(cell_num(1)*d(2)/d(1))  + 1, 1000)
            else if (dim == 3) then
                cell_num(1) = min(int(((ntotal)*d(1)*d(1) &
                                /(d(2)*d(3)*nppg))**(1._8/3)) + 1, 1000)
                cell_num(2) = min(int(cell_num(1)*d(2)/d(1))  + 1, 1000)
                cell_num(3) = min(int(cell_num(1)*d(3)/d(1))  + 1, 1000)
            end if
            !!! smoothing length measured in grid cells:
            ghsmlx(:) = int(dble(cell_num(:))*hsml/d(:)) + 1
        end associate

        !!! Initialize grid
        select case(dim)
        case (1)
            allocate( grid(cell_num(1),          1,           1 ), source=0 )
        case (2)
            allocate( grid(cell_num(1), cell_num(2),          1 ), source=0 )
        case (3)
            allocate( grid(cell_num(1), cell_num(2), cell_num(3)), source=0 )
        end select

    end subroutine initialize_grid

    subroutine grid_geometry(index, coor, cell_num, grid_max_coor, grid_min_coor, cell)
        use ctrl_dict, only: i_time_step
        integer, intent(in) :: index            !! Particle index
        real(8), intent(in) :: coor(:)          !! Particle coordinates
        integer, intent(in) :: cell_num(:)      !! Cell number
        real(8), intent(in) :: grid_max_coor(:), &
                               grid_min_coor(:) !! Maximum and minium grid cell coordinates
        integer, intent(inout) :: cell(3)     !! Cell coordinates

        integer :: d

        cell(:) = 1

        do d = 1, dim
            if ((coor(d) > grid_max_coor(d)) .or. (coor(d) < grid_min_coor(d))) then
                write(*,*) ' >> error: particle out of range'
                write(*,*) "   Time Step:"//to_string(i_time_step)
                write(*,*) '   particle position: coor('//to_string(index)//', ' &
                            //to_string(d)//') = '//to_string(coor(d))
                write(*,*) '   range: [xmin,xmax]('//to_string(d)//') = ['// &
                            to_string(grid_min_coor(d))//', '//to_string(grid_max_coor(d))//']'
                error stop
            else
                cell(d) = int(dble(cell_num(d))*(coor(d) - grid_min_coor(d))/(grid_max_coor(d)-grid_min_coor(d)) + 1.0_8)
                if (cell(d) > cell_num(d)) then
                    write(*,*) 'cell(d) is greater than cell_num(d)'
                    cell(d) = cell_num(d)
                end if
            end if
        end do

    end subroutine grid_geometry

    subroutine tree_search(P)
        use tree_m
        use link_list_m
        use geometry_m
        type(Particle), intent(inout) :: P(:)
        integer :: ntotal, kpair
        real(8) :: this_w, this_dwdx(dim)
        real(8) :: mhsml
        class(geometry_t), allocatable :: domain, range
        type(tree_t)      :: tree
        type(link_list_t) :: found
        class(*), allocatable :: j
        real(8) :: scale, min(dim), max(dim), length(dim)
        real(8) :: dx(dim), dr, r
        integer :: scale_k
        logical :: flag
        integer i, k, d

        found%length = 0

        select case (skf)
        case (1)
            scale_k = 2
        case (2, 3)
            scale_k = 3
        case default
            scale_k = 1
        end select

        ntotal = size(P)
        kpair = size(P(1)%neighborList)

        scale = 1.2
        min = huge(0._8)
        max = -huge(0._8)
        do i = 1, ntotal
            do d = 1, dim
                if ( P(i)%x(d) > max(d) ) then
                    max(d) = P(i)%x(d)
                end if
                if ( P(i)%x(d) < min(d) ) then
                    min(d) = P(i)%x(d)
                end if
            end do
        end do
        length = scale * (max - min)
        select case (dim)
        case (1)
            domain = line_t((min+max)/2, length(1), 0)
        case (2)
            domain = rectangle_t((min+max)/2, length, 0)
        case (3)
            domain = cuboid_t((min+max)/2, length, 0)
        end select

        tree%domain   = domain
        tree%capacity = 4
        allocate(tree%points(0))
        do i = 1, ntotal
            call tree%insert(point_t(P(i)%x(:), i), flag)
            if ( .not. flag ) then
                write(err, *) i
                error stop "Tree insert"
            end if
        end do

        deallocate(domain)

        ! !$OMP PARALLEL DO PRIVATE(i, j, k, d, range, found, mhsml, dx, dr, r, this_w, this_dwdx) &
        ! !$OMP SHARED(P) &
        ! !$OMP SCHEDULE(dynamic, chunkSize)
        do i = 1, ntotal - 1
            select case (dim)
            case (1)
                range = line_t(P(i)%x(:), scale_k*P(i)%SmoothingLength*2, 0)
            case (2)
                range = circle_t(P(i)%x(:), scale_k*P(i)%SmoothingLength, 0)
            case (3)
                range = sphere_t(P(i)%x(:), scale_k*P(i)%SmoothingLength, 0)
            end select
            ! !$OMP CRITICAL
            call tree%query(range, found)
            ! !$OMP END CRITICAL
            do k = 1, found%length
                call found%fetch(j)
                select type (j)
                type is (integer)
                    if ( j > i ) then
                        mhsml = 0.5_8 * (P(i)%SmoothingLength + P(j)%SmoothingLength)
                        dx(1) = P(i)%x(1) - P(j)%x(1)
                        dr = dx(1) * dx(1)
                        do d = 2, dim
                            dx(d) = P(i)%x(d) - P(j)%x(d)
                            dr = dr + dx(d)*dx(d)
                        end do
                        r = sqrt(dr)
                        ! !$OMP CRITICAL
                        !!! Neighbor particle ant tottal neighbor number for each particle
                        P(i)%neighborNum = P(i)%neighborNum + 1
                        P(j)%neighborNum = P(j)%neighborNum + 1
#if CHECK_NEIGHBOR_NUM
                        if ( P(i)%neighborNum > kpair .or. P(j)%neighborNum > kpair ) then
                            write(*,*) i, j, P([i,j])%neighborNum, kpair
                            error stop "Too many neighbors"
                        end if
#endif
                        P(i)%neighborList(P(i)%neighborNum) = j
                        P(j)%neighborList(P(j)%neighborNum) = i
                        !!! Kernel and derivations of kernel
                        call kernel(r, dx, mhsml, this_w, this_dwdx)
                        P(i)%w(P(i)%neighborNum) = this_w
                        P(j)%w(P(j)%neighborNum) = this_w
                        P(i)%dwdx(:, P(i)%neighborNum) =   this_dwdx
                        P(j)%dwdx(:, P(j)%neighborNum) = - this_dwdx
                        ! !$OMP END CRITICAL
                    end if
                class default
                    error stop "Tree query"
                end select
                ! deallocate(j)
            end do
            ! deallocate(range)
        end do
        ! !$OMP END PARALLEL DO

        call tree%clean()

    end subroutine tree_search


    subroutine print_statistics(Particles)
        type(Particle), intent(in) :: Particles(:)
        integer :: ntotal
        integer :: sumiac, maxiac, miniac, noiac, &
                   maxp(1), minp(1)
        integer i

        sumiac =  0
        maxiac =  0
        miniac =  huge(0)
        noiac  =  0
        maxp   = -1
        minp   = -1

        ntotal = size(Particles)
        do i = 1, ntotal
            sumiac = sumiac + Particles(i)%neighborNum
            if (Particles(i)%neighborNum > maxiac) then
                maxiac = Particles(i)%neighborNum
                maxp = i
            end if
            if (Particles(i)%neighborNum < miniac) then
                miniac = Particles(i)%neighborNum
                minp = i
            end if
            if (Particles(i)%neighborNum == 0) noiac = noiac + 1
        end do

        write(*,*) ">> Statistics: interactions per particles:"
        write(*,*) "   Particle "//to_string(maxp(1))// &
                   " has maximal interactions: "//to_string(maxiac)
        write(*,*) "   Particle "//to_string(minp(1))// &
                   " has minimal interactions: "//to_string(miniac)
        write(*, "(A,G0)") "    Average neighborList per particle: ", &
                           real(sumiac)/(ntotal)
        write(*,*) "   Particle with no interacitons: "//to_string(noiac)



    end subroutine print_statistics

end module nnps_m