module tree_m
    use geometry_m
    use link_list_m
    implicit none
    private

    type, public :: tree_t
        integer :: capacity
        logical :: splited = .false.
        type(point_t),     allocatable :: points(:)
        class(geometry_t), allocatable :: domain
        type(tree_t), pointer     :: parent
        type(tree_t), allocatable :: children(:)
        integer :: rank
    contains
        procedure :: insert => subtree_insert
        procedure :: split  => subtree_split
        procedure :: query  => subtree_query
        procedure :: clean  => subtree_clean
    end type tree_t
    
contains
    recursive pure subroutine subtree_insert(this, point, flag)
        class(tree_t), intent(inout) :: this
        type(point_t), intent(in)    :: point
        logical,       intent(inout) :: flag
        integer :: i

        flag = .false.
        if (this%domain%contain(point)) then
            if ( size(this%points) < this%capacity ) then
                this%points = [this%points, point]
                flag = .true.
            else
                if (.not. this%splited) call this%split()
                do i = 1, size(this%children)
                    call this%children(i)%insert(point, flag)
                    if (flag) return
                end do
            end if
        end if

    end subroutine subtree_insert

    pure subroutine subtree_split(this)
        class(tree_t), intent(inout) :: this
        real(8) :: scale
        ! integer i

        scale = 0.5 + 0.001

        associate(domain => this%domain)
        select type (domain)
        type is (line_t)
            allocate(this%children(2))
            this%children(:)%capacity = this%capacity
            select case (size(domain%center))
            case (1)
                this%children(1)%domain = line_t(domain%center - domain%length/4, &
                                                 scale*domain%length,             &
                                                 domain%index + 1)
                this%children(2)%domain = line_t(domain%center + domain%length/4, &
                                                 scale*domain%length,             &
                                                 domain%index + 1)
                allocate(this%children(1)%points(0))
                allocate(this%children(2)%points(0))
            case (2)
                this%children(1)%domain = line_t([domain%center(1) - domain%length/4], &
                                                  scale*domain%length,                 &
                                                  domain%index + 1)
                this%children(2)%domain = line_t([domain%center(1) + domain%length/4], &
                                                  scale*domain%length,                 &
                                                  domain%index + 1)
                allocate(this%children(1)%points(0))
                allocate(this%children(2)%points(0))
            case (3)
                this%children(1)%domain = line_t([domain%center(1) - domain%length/4], &
                                                  scale*domain%length,                 &
                                                  domain%index + 1)
                this%children(2)%domain = line_t([domain%center(1) + domain%length/4], &
                                                  scale*domain%length,                 &
                                                  domain%index + 1)
                allocate(this%children(1)%points(0))
                allocate(this%children(2)%points(0))
            end select
        type is (rectangle_t)
            allocate(this%children(4))
            this%children(:)%capacity = this%capacity
            select case (size(domain%center))
            case (2)
                this%children(1)%domain = rectangle_t([domain%center(1) - domain%length(1)/4,  &
                                                       domain%center(2) - domain%length(2)/4], &
                                                       scale*domain%length,                    &
                                                       domain%index + 1)
                this%children(2)%domain = rectangle_t([domain%center(1) + domain%length(1)/4,  &
                                                       domain%center(2) - domain%length(2)/4], &
                                                       scale*domain%length,                    &
                                                       domain%index + 1)
                this%children(3)%domain = rectangle_t([domain%center(1) - domain%length(1)/4,  &
                                                       domain%center(2) + domain%length(2)/4], &
                                                       scale*domain%length,                    &
                                                       domain%index + 1)
                this%children(4)%domain = rectangle_t([domain%center(1) + domain%length(1)/4,  &
                                                       domain%center(2) + domain%length(2)/4], &
                                                       scale*domain%length,                    &
                                                       domain%index + 1)
                allocate(this%children(1)%points(0))
                allocate(this%children(2)%points(0))
                allocate(this%children(3)%points(0))
                allocate(this%children(4)%points(0))
            case (3)
                this%children(1)%domain = rectangle_t([domain%center(1) - domain%length(1)/4,  &
                                                       domain%center(2) - domain%length(2)/4], &
                                                       scale*domain%length,                    &
                                                       domain%index + 1)
                this%children(2)%domain = rectangle_t([domain%center(1) + domain%length(1)/4,  &
                                                       domain%center(2) - domain%length(2)/4], &
                                                       scale*domain%length,                    &
                                                       domain%index + 1)
                this%children(3)%domain = rectangle_t([domain%center(1) - domain%length(1)/4,  &
                                                       domain%center(2) + domain%length(2)/4], &
                                                       scale*domain%length,                    &
                                                       domain%index + 1)
                this%children(4)%domain = rectangle_t([domain%center(1) + domain%length(1)/4,  &
                                                       domain%center(2) + domain%length(2)/4], &
                                                       scale*domain%length,                    &
                                                       domain%index + 1)
                allocate(this%children(1)%points(0))
                allocate(this%children(2)%points(0))
                allocate(this%children(3)%points(0))
                allocate(this%children(4)%points(0))
            end select
        type is (cuboid_t)
            allocate(this%children(8))
            this%children(:)%capacity = this%capacity
            this%children(1)%domain = cuboid_t([domain%center(1) - domain%length(1)/4,  &
                                                domain%center(2) - domain%length(2)/4,  &
                                                domain%center(3) - domain%length(3)/4], &
                                                scale*domain%length,                    &
                                                domain%index + 1)
            this%children(2)%domain = cuboid_t([domain%center(1) + domain%length(1)/4,  &
                                                domain%center(2) - domain%length(2)/4,  &
                                                domain%center(3) - domain%length(3)/4], &
                                                scale*domain%length,                    &
                                                domain%index + 1)
            this%children(3)%domain = cuboid_t([domain%center(1) - domain%length(1)/4,  &
                                                domain%center(2) + domain%length(2)/4,  &
                                                domain%center(3) - domain%length(3)/4], &
                                                scale*domain%length,                    &
                                                domain%index + 1)
            this%children(4)%domain = cuboid_t([domain%center(1) + domain%length(1)/4,  &
                                                domain%center(2) + domain%length(2)/4,  &
                                                domain%center(3) - domain%length(3)/4], &
                                                scale*domain%length,                    &
                                                domain%index + 1)
            this%children(5)%domain = cuboid_t([domain%center(1) - domain%length(1)/4,  &
                                                domain%center(2) - domain%length(2)/4,  &
                                                domain%center(3) + domain%length(3)/4], &
                                                scale*domain%length,                    &
                                                domain%index + 1)
            this%children(6)%domain = cuboid_t([domain%center(1) + domain%length(1)/4,  &
                                                domain%center(2) - domain%length(2)/4,  &
                                                domain%center(3) + domain%length(3)/4], &
                                                scale*domain%length,                    &
                                                domain%index + 1)
            this%children(7)%domain = cuboid_t([domain%center(1) - domain%length(1)/4,  &
                                                domain%center(2) + domain%length(2)/4,  &
                                                domain%center(3) + domain%length(3)/4], &
                                                scale*domain%length,                    &
                                                domain%index + 1)
            this%children(8)%domain = cuboid_t([domain%center(1) + domain%length(1)/4,  &
                                                domain%center(2) + domain%length(2)/4,  &
                                                domain%center(3) + domain%length(3)/4], &
                                                scale*domain%length,                    &
                                                domain%index + 1)
                allocate(this%children(1)%points(0))
                allocate(this%children(2)%points(0))
                allocate(this%children(3)%points(0))
                allocate(this%children(4)%points(0))
                allocate(this%children(5)%points(0))
                allocate(this%children(6)%points(0))
                allocate(this%children(7)%points(0))
                allocate(this%children(8)%points(0))
        end select
        end associate
        this%splited = .true.
    end subroutine subtree_split

    recursive pure subroutine subtree_query(this, range, found)
        class(tree_t), intent(inout) :: this
        class(geometry_t),  intent(in)    :: range
        class(link_list_t), intent(inout) :: found
        integer i

        if (range%intersect(this%domain)) then
            do i = 1, size(this%points)
                if (range%contain(this%points(i))) then
                    call found%add(this%points(i)%index)
                end if
            end do
            if (this%splited) then
                do i = 1, size(this%children)
                    call this%children(i)%query(range, found)
                end do
            end if
        end if
        
    end subroutine subtree_query

    recursive pure subroutine subtree_clean(this)
        class(tree_t), intent(inout) :: this
        integer i

        ! call this%domain%clean()
        deallocate(this%domain)
        if ( allocated(this%points) .and. size(this%points) > 0 ) then
            deallocate(this%points)
        end if
        if ( this%splited ) then
            do i = 1, size(this%children)
                call this%children(i)%clean()
            end do
        end if

    end subroutine
    
end module tree_m