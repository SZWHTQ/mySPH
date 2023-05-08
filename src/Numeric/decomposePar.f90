module decompose_m
    implicit none

contains

    subroutine factor2(input, output, count)
        implicit none
        integer, intent(in)  :: input
        integer, intent(out) :: output(2)
        integer, intent(out), optional :: count

        if (present(count)) count = 1
        output(2) = int(sqrt(dble(input)))
        output(1) = int(dble(input) / output(2))
        do while( abs(output(1) * output(2) - input) /= 0 )
            output(2) = output(2) - 1
            output(1) = int(dble(input) / output(2))
            if (present(count)) count = count + 1
        end do

    end subroutine factor2

    subroutine factor3(input, output, count)
        implicit none
        integer, intent(in)  :: input
        integer, intent(out) :: output(3)
        integer, intent(out), optional :: count
        integer :: buffer

        if (present(count)) count = 1
        output(3) = int(dble(input)**(1./3))
        output(2) = int(dble(input) / output(3))
        do while( abs(output(3) * output(2) - input) /= 0 )
            output(3) = output(3) - 1
            output(2) = int(dble(input) / output(3))
            if (present(count)) count = count + 1
        end do
        call factor2(output(2)*1, output(1:2), buffer)
        if (present(count)) count = count + buffer

    end subroutine factor3

    subroutine decompose(Type, x, extend, nsub, sub_ntotal, sub_ndummy, sub_index)
        use geometry_m
        integer, intent(in)  :: Type(:)
        real(8), intent(in)  :: x(:, :)
        real(8), intent(in)  :: extend
        integer, intent(in)  :: nsub
        integer, intent(out) :: sub_ntotal(:)
        integer, intent(out) :: sub_ndummy(:)
        integer, intent(out) :: sub_index(:, :)
        integer, allocatable :: sub(:)
        class(geometry_t), allocatable :: grid(:, :, :, :)
        save grid
        integer :: dim
        ! integer i, j, k

        dim = size(x, 1)
        allocate(sub(dim))

        if ( .not. allocated(grid) )  then
            call initialize_grid(dim, nsub, sub, grid)
            call equally_decompose(x, extend, grid)
            ! write(*,*) "Initialize Grid"
        end if

        call distribute_point(Type, x, grid, sub_ntotal, sub_ndummy, sub_index)

        ! write(*, "(A)") " >> Statistics: OpenMP Grid Info: "
        ! write(*, "(3X, A15, 3(A17), 3(A8))") "Center-X", "Center-Y", "Length-X", "Length-Y", &
        !                             "Index", "Ntotal", "Ndummy"
        ! select type (grid)
        ! type is (line_t)
        !     do i = 1, size(grid, 1)
        !         write(*, "(3X, 2(E15.7, 2X), 3(6X, I0))") grid(i, 1, 1, 1), sub_ntotal(i), sub_ndummy(i)
        !     end do
        ! type is (rectangle_t)
        !     do i = 1, size(grid, 1)
        !         do j = 1, size(grid, 2)
        !             write(*, "(3X, 4(E15.7, 2X), 3(6X, I0))") grid(i, j, 1, 1), sub_ntotal(i), sub_ndummy(i)
        !         end do
        !     end do
        ! type is (cuboid_t)
        !     do i = 1, size(grid, 1)
        !         do j = 1, size(grid, 2)
        !             do k = 1, size(grid, 3)
        !                 write(*, "(3X, 6(E15.7, 2X), 3(6X, I0))") grid(i, j, k, 1), sub_ntotal(i), sub_ndummy(i)
        !             end do
        !         end do
        !     end do
        ! end select

        deallocate(grid)

    end subroutine decompose

    subroutine initialize_grid(dim, nsub, sub, grid)
        use geometry_m
        integer, intent(in)  :: dim
        integer, intent(in)  :: nsub
        integer, intent(out) :: sub(:)
        class(geometry_t), allocatable, intent(out) :: grid(:, :, :, :)

        select case (dim)
        case(1)
            allocate(grid(nsub,   1,      1, 2), &
                     source = line_t([0._8], 0._8, 0))
        case (2)
            call factor2(nsub, sub)
            allocate(grid(sub(1), sub(2), 1, 2), &
                     source = rectangle_t([0._8, 0._8], 0._8, 0))
        case (3)
            call factor3(nsub, sub)
            allocate(grid(sub(1), sub(2), sub(3), 2), &
                     source = cuboid_t([0._8, 0._8, 0._8], 0._8, 0))
        end select

    end subroutine initialize_grid

    subroutine equally_decompose(x, extend, grid)
        use geometry_m
        real(8), intent(in)  :: x(:, :)
        real(8), intent(in)  :: extend
        class(geometry_t), intent(out) :: grid(:, :, :, :)
        integer :: dim, num
        real(8), allocatable :: min_coor(:), max_coor(:)
        real(8), allocatable :: length(:)
        integer i, j, k, d

        dim = size(x, 1)
        allocate(min_coor(dim), max_coor(dim))
        allocate(length(dim))

        do d = 1, dim
            num = size(grid, d)
            min_coor(d) = minval(x(d, :))
            max_coor(d) = maxval(x(d, :))
            length(d) = (max_coor(d) - min_coor(d)) / num
        end do

        select type (grid)
        type is (line_t)
            do i = 1, size(grid, 1)
                grid(i, 1, 1, 1) = line_t(min_coor+length*(real(i)-0.5), &
                                          length(1),                     &
                                          i)
                grid(i, 1, 1, 2) = line_t(min_coor+length*(real(i)-0.5), &
                                          length(1)+extend,              &
                                          i)
            end do
        type is (rectangle_t)
            do i = 1, size(grid, 1)
                do j = 1, size(grid, 2)
                    grid(i, j, 1, 1) = rectangle_t(min_coor+length*(real([i, j])-0.5), &
                                                   length,                             &
                                                   (i-1)*size(grid, 2) + j)
                    grid(i, j, 1, 2) = rectangle_t(min_coor+length*(real([i, j])-0.5), &
                                                   length+extend,                      &
                                                   (i-1)*size(grid, 2) + j)
                end do
            end do
        type is (cuboid_t)
            do i = 1, size(grid, 1)
                do j = 1, size(grid, 2)
                    do k = 1, size(grid, 3)
                        grid(i, j, k, 1) = cuboid_t(min_coor+length*(real([i, j, k])-0.5), &
                                                    length,                                &
                                                    (i-1)*size(grid, 2)*size(grid, 3) + (j-1)*size(grid, 3) + k)
                        grid(i, j, k, 2) = cuboid_t(min_coor+length*(real([i, j, k])-0.5), &
                                                    length+extend,                         &
                                                    (i-1)*size(grid, 2)*size(grid, 3) + (j-1)*size(grid, 3) + k)
                    end do
                end do
            end do
        end select

    end subroutine equally_decompose

    subroutine distribute_point(Type, x, grid, sub_ntotal, sub_ndummy, sub_index)
        !$ use omp_lib
        use geometry_m
        integer, intent(in)  :: Type(:)
        real(8), intent(in)  :: x(:, :)
        class(geometry_t), intent(in)  :: grid(:, :, :, :)
        integer, intent(out) :: sub_ntotal(:)
        integer, intent(out) :: sub_ndummy(:)
        integer, intent(out) :: sub_index(:, :)
        logical :: distributed
        integer :: cell_id
        integer n, i, j, k

        sub_ntotal = 0
        sub_ndummy = 0
        sub_index  = 0
        distributed = .false.

        select type (grid)
        type is (line_t)
            !!$OMP PARALLEL DO PRIVATE(n, i, cell_id, distributed)
            do n = 1, size(x, 2)
                distributed = .false.
                do i = 1, size(grid, 1)
                    if ( grid(i, 1, 1, 1)%contain(point_t(x(:,n), n)) ) then
                        cell_id = grid(i, 1, 1, 1)%index
                        if ( Type(n) > 0 ) then
                            sub_ntotal(cell_id) = sub_ntotal(cell_id) + 1
                            sub_index(cell_id,sub_ntotal(cell_id)) = n
                        end if
                        distributed = .true.
                    end if
                    if ( distributed ) exit
                end do
            end do
            !!$OMP END PARALLEL DO
        type is (rectangle_t)
            !!$OMP PARALLEL DO PRIVATE(n, i, j, cell_id, distributed)
            do n = 1, size(x, 2)
                distributed = .false.
                do i = 1, size(grid, 1)
                    do j = 1, size(grid, 2)
                        if ( grid(i, j, 1, 1)%contain(point_t(x(:,n), n)) ) then
                            cell_id = grid(i, j, 1, 1)%index
                            if ( Type(n) > 0 ) then
                                sub_ntotal(cell_id) = sub_ntotal(cell_id) + 1
                                sub_index(cell_id,sub_ntotal(cell_id)) = n
                            end if
                            distributed = .true.
                        end if
                        if ( distributed ) exit
                    end do
                    if ( distributed ) exit
                end do
            end do
            !!$OMP END PARALLEL DO
        type is (cuboid_t)
            !!$OMP PARALLEL DO PRIVATE(n, i, j, k, cell_id, distributed)
            do n = 1, size(x, 2)
                distributed = .false.
                do i = 1, size(grid, 1)
                    do j = 1, size(grid, 2)
                        do k = 1, size(grid, 3)
                            if ( grid(i, j, k, 1)%contain(point_t(x(:,n), n)) ) then
                                cell_id = grid(i, j, k, 1)%index
                                if ( Type(n) > 0 ) then
                                    sub_ntotal(cell_id) = sub_ntotal(cell_id) + 1
                                    sub_index(cell_id,sub_ntotal(cell_id)) = n
                                end if
                                distributed = .true.
                            end if
                            if ( distributed ) exit
                        end do
                        if ( distributed ) exit
                    end do
                    if ( distributed ) exit
                end do
            end do
            !!$OMP END PARALLEL DO
        end select

        select type (grid)
        type is (line_t)
            !!$OMP PARALLEL DO PRIVATE(n, i, cell_id, distributed)
            do n = 1, size(x, 2)
                distributed = .false.
                do i = 1, size(grid, 1)
                    if ( grid(i, 1, 1, 1)%contain(point_t(x(:,n), n)) ) then
                        if ( Type(n) < 0 ) then
                            cell_id = grid(i, 1, 1, 1)%index
                            sub_ndummy(cell_id) = sub_ndummy(cell_id) + 1
                            sub_index(cell_id,sub_ntotal(cell_id)+sub_ndummy(cell_id)) = n
                        end if
                        distributed = .true.
                    end if
                    if ( distributed ) exit
                end do
            end do
            !!$OMP END PARALLEL DO
        type is (rectangle_t)
            !!$OMP PARALLEL DO PRIVATE(n, i, j, cell_id, distributed)
            do n = 1, size(x, 2)
                distributed = .false.
                do i = 1, size(grid, 1)
                    do j = 1, size(grid, 2)
                        if ( grid(i, j, 1, 1)%contain(point_t(x(:,n), n)) ) then
                            if ( Type(n) < 0 ) then
                                cell_id = grid(i, j, 1, 1)%index
                                sub_ndummy(cell_id) = sub_ndummy(cell_id) + 1
                                sub_index(cell_id,sub_ntotal(cell_id)+sub_ndummy(cell_id)) = n
                            end if
                            distributed = .true.
                        end if
                        if ( distributed ) exit
                    end do
                    if ( distributed ) exit
                end do
            end do
            !!$OMP END PARALLEL DO
        type is (cuboid_t)
            !!$OMP PARALLEL DO PRIVATE(n, i, j, k, cell_id, distributed)
            do n = 1, size(x, 2)
                distributed = .false.
                do i = 1, size(grid, 1)
                    do j = 1, size(grid, 2)
                        do k = 1, size(grid, 3)
                                if ( grid(i, j, k, 1)%contain(point_t(x(:,n), n)) ) then
                                    if ( Type(n) < 0 ) then
                                        cell_id = grid(i, j, k, 1)%index
                                        sub_ndummy(cell_id) = sub_ndummy(cell_id) + 1
                                        sub_index(cell_id,sub_ntotal(cell_id)+sub_ndummy(cell_id)) = n
                                    end if
                                    distributed = .true.
                                end if
                            if ( distributed ) exit
                        end do
                        if ( distributed ) exit
                    end do
                    if ( distributed ) exit
                end do
            end do
            !!$OMP END PARALLEL DO
        end select

        select type (grid)
        type is (line_t)
            !!$OMP PARALLEL DO PRIVATE(n, i, cell_id, distributed)
            do n = 1, size(x, 2)
                do i = 1, size(grid, 1)
                    if ( grid(i, 1, 1, 2)%contain(point_t(x(:,n), n)) ) then
                        if ( .not. grid(i, 1, 1, 1)%contain(point_t(x(:,n), n)) ) then
                            cell_id = grid(i, 1, 1, 1)%index
                            sub_ndummy(cell_id) = sub_ndummy(cell_id) + 1
                            sub_index(cell_id,sub_ntotal(cell_id)+sub_ndummy(cell_id)) = n
                        end if
                    end if
                end do
            end do
            !!$OMP END PARALLEL DO
        type is (rectangle_t)
            !!$OMP PARALLEL DO PRIVATE(n, i, j, cell_id, distributed)
            do n = 1, size(x, 2)
                do i = 1, size(grid, 1)
                    do j = 1, size(grid, 2)
                        if ( grid(i, j, 1, 2)%contain(point_t(x(:,n), n)) ) then
                            if ( .not. grid(i, j, 1, 1)%contain(point_t(x(:,n), n)) ) then
                                cell_id = grid(i, j, 1, 1)%index
                                sub_ndummy(cell_id) = sub_ndummy(cell_id) + 1
                                sub_index(cell_id,sub_ntotal(cell_id)+sub_ndummy(cell_id)) = n
                            end if
                        end if
                    end do
                end do
            end do
            !!$OMP END PARALLEL DO
        type is (cuboid_t)
            !!$OMP PARALLEL DO PRIVATE(n, i, j, k, cell_id, distributed)
            do n = 1, size(x, 2)
                do i = 1, size(grid, 1)
                    do j = 1, size(grid, 2)
                        do k = 1, size(grid, 3)
                            if ( grid(i, j, k, 2)%contain(point_t(x(:,n), n)) ) then
                                if ( .not. grid(i, j, k, 1)%contain(point_t(x(:,n), n)) ) then
                                    cell_id = grid(i, j, k, 1)%index
                                    sub_ndummy(cell_id) = sub_ndummy(cell_id) + 1
                                    sub_index(cell_id,sub_ntotal(cell_id)+sub_ndummy(cell_id)) = n
                                end if
                            end if
                        end do
                    end do
                end do
            end do
            !!$OMP END PARALLEL DO
        end select

    end subroutine distribute_point

end module decompose_m