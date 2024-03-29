module back_ground_grid_m
    use ctrl_dict, only: Field
    implicit none
    
contains
    pure subroutine initialize_grid(ntotal, hsml, grid, cell_num, ghsmlx, gridMaxCoor, gridMinCoor)
        integer, intent(in) :: ntotal                                !! Number of particles
        real(8), intent(in) :: hsml                                  !! Smoothing length
        integer, intent(inout), allocatable :: grid(:, :, :)         !! Grid
        integer, intent(inout) :: cell_num(:)                        !! Number of grid cells
        integer, intent(inout) :: ghsmlx(:)                          !! Smoothing length
        real(8), intent(inout) :: gridMaxCoor(:), gridMinCoor(:) !! Maximum and minium grid cell coordinates
        integer, parameter     :: nppg = 3  !! averaged number of particles per grid cell

        !!! initialize parameters: maximum number of grid cells

        !!! range of sorting grid

        associate( d=>(gridMaxCoor(:) - gridMinCoor(:)) )
            !!! number of grid cells in x-, y- and z-direction:
            if (Field%Dim == 1) then
                cell_num(1) = min((ntotal)/nppg + 1, 1000)
            else if (Field%Dim == 2) then
                cell_num(1) = min(int(((ntotal) * d(1) &
                                / (d(2)*nppg))**(1._8/2) )   + 1, 1000)
                cell_num(2) = min(int(cell_num(1)*d(2)/d(1)) + 1, 1000)
            else if (Field%Dim == 3) then
                cell_num(1) = min(int(((ntotal) * d(1) * d(1) &
                                / (d(2) * d(3) * nppg))**(1._8/3))+ 1, 1000)
                cell_num(2) = min(int(cell_num(1) * d(2) / d(1))  + 1, 1000)
                cell_num(3) = min(int(cell_num(1) * d(3) / d(1))  + 1, 1000)
            end if
            !!! smoothing length measured in grid cells:
            ghsmlx(:) = int(dble(cell_num(:))*hsml/d(:)) + 1
        end associate

        !!! Initialize grid
        select case(Field%Dim)
        case (1)
            allocate( grid(cell_num(1),          1,           1 ), source=0 )
        case (2)
            allocate( grid(cell_num(1), cell_num(2),          1 ), source=0 )
        case (3)
            allocate( grid(cell_num(1), cell_num(2), cell_num(3)), source=0 )
        end select

    end subroutine initialize_grid

    subroutine grid_geometry(index, coor, cell_num, gridMaxCoor, gridMinCoor, cell)
        use ctrl_dict, only: Config
        use tools_m,   only: to_string
        integer, intent(in) :: index            !! Particle index
        real(8), intent(in) :: coor(:)          !! Particle coordinates
        integer, intent(in) :: cell_num(:)      !! Cell number
        real(8), intent(in) :: gridMaxCoor(:), &
                            gridMinCoor(:) !! Maximum and minium grid cell coordinates
        integer, intent(inout) :: cell(3)     !! Cell coordinates

        integer :: d

        cell(:) = 1

        do d = 1, Field%Dim
            if ((coor(d) > gridMaxCoor(d)) .or. (coor(d) < gridMinCoor(d))) then
                write(*,*) ' >> error: particle out of range'
                write(*,*) "   Time Step:"//to_string(Config%i_time_step)
                write(*,*) '   particle position: coor('//to_string(index)//', ' &
                            //to_string(d)//') = '//to_string(coor(d))
                write(*,*) '   range: [xmin,xmax]('//to_string(d)//') = ['// &
                            to_string(gridMinCoor(d))//', '//to_string(gridMaxCoor(d))//']'
                error stop
            else
                cell(d) = int(dble(cell_num(d))*(coor(d) - gridMinCoor(d))/(gridMaxCoor(d)-gridMinCoor(d)) + 1.0_8)
                if (cell(d) > cell_num(d)) then
                    write(*,*) 'cell(d) is greater than cell_num(d)'
                    cell(d) = cell_num(d)
                end if
            end if
        end do

    end subroutine grid_geometry

end module back_ground_grid_m