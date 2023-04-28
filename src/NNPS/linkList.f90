module link_list_m
    implicit none
    private

    type, public :: link_node_t
        class(*), allocatable :: item
        ! type(link_node_t), pointer :: prev => null()
        type(link_node_t), pointer :: next => null()
    end type link_node_t

    type, public :: link_list_t
        type(link_node_t), pointer :: head => null()
        type(link_node_t), pointer :: tail => null()
        integer :: length
    contains
        procedure :: add   => link_list_add
        procedure :: fetch => link_list_fetch
    end type link_list_t

contains
    ! pure function init_node(item) result(node)
    !         class(*), intent(in) :: item
    !         type(link_node_t) :: node
    !         allocate(node%item, source=item)
    ! end function init_node

    pure subroutine link_list_add(this, item)
        class(link_list_t), intent(inout) :: this
        class(*), intent(in) :: item
        type(link_node_t) :: node

        if ( associated(this%tail) ) then
            allocate(node%item,      source=item)
            allocate(this%tail%next, source=node)
            this%tail => this%tail%next
        else
            allocate(node%item, source=item)
            allocate(this%head, source=node)
            this%tail => this%head
        end if
        this%length = this%length + 1

    end subroutine link_list_add

    pure subroutine link_list_fetch(this, item)
        class(link_list_t), intent(inout) :: this
        class(*), intent(inout), allocatable, optional :: item
        type(link_node_t), pointer :: node

        if ( associated(this%head) ) then
            if ( present(item) ) then
                call move_alloc(this%head%item, item)
            else
                deallocate(this%head%item)
            end if
            node      => this%head
            this%head => this%head%next
            this%length = this%length - 1
            nullify(node%next)
            deallocate(node)
            if ( this%length == 0 ) then
                nullify(this%head, this%tail)
            end if
        end if

    end subroutine link_list_fetch

end module link_list_m