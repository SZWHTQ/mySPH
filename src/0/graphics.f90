module geometry_m
    use, intrinsic :: iso_fortran_env, only: err => error_unit
    implicit none
    private

    type, public, abstract :: geometry_t
    contains
        procedure(graphics_contain),   deferred :: contain
        procedure(graphics_intersect), deferred :: intersect
        ! procedure(graphics_clean),     deferred :: clean
    end type geometry_t

    abstract interface
        logical elemental function graphics_contain(this, that)
            import geometry_t
            class(geometry_t), intent(in) :: this, that
        end function graphics_contain
        logical elemental function graphics_intersect(this, that)
            import geometry_t
            class(geometry_t), intent(in) :: this, that
        end function graphics_intersect
        ! pure subroutine graphics_clean(this)
        !     import geometry_t
        !     class(geometry_t), intent(inout) :: this
        ! end subroutine graphics_clean
    end interface

    type, public, extends(geometry_t) :: point_t
        real(8), allocatable :: center(:)
        integer :: index
    contains
        procedure :: contain   => point_contain
        procedure :: intersect => point_intersect
        ! procedure :: clean     => point_clean
    end type point_t

    type, public, extends(geometry_t) :: line_t
        real(8) :: center(1)
        real(8) :: length
        integer :: index
    contains
        procedure :: contain   => line_contain
        procedure :: intersect => line_intersect
        ! procedure :: clean     => line_clean
    end type line_t

    type, public, extends(geometry_t) :: rectangle_t
        real(8) :: center(2)
        real(8) :: length(2)
        integer :: index
    contains
        procedure :: contain   => rectangle_contain
        procedure :: intersect => rectangle_intersect
        ! procedure :: clean     => rectangle_clean
    end type rectangle_t

    type, public, extends(geometry_t) :: circle_t
        real(8) :: center(2)
        real(8) :: radius
        integer :: index
    contains
        procedure :: contain   => circle_contain
        procedure :: intersect => circle_intersect
        ! procedure :: clean     => circle_clean
    end type circle_t

    type, public, extends(geometry_t) :: cuboid_t
        real(8) :: center(3)
        real(8) :: length(3)
        integer :: index
    contains
        procedure :: contain   => cuboid_contain
        procedure :: intersect => cuboid_intersect
        ! procedure :: clean     => cuboid_clean
    end type cuboid_t

    type, public, extends(geometry_t) :: sphere_t
        real(8) :: center(3)
        real(8) :: radius
        integer :: index
    contains
        procedure :: contain   => sphere_contain
        procedure :: intersect => sphere_intersect
        ! procedure :: clean     => sphere_clean
    end type sphere_t

contains
    logical elemental function point_contain(this, that) result(result)
        class(point_t), intent(in)    :: this
        class(geometry_t), intent(in) :: that
        select type(that)
        type is (point_t)
            result = all(this%center == that%center)
        class default
            result = .false.
        end select
    end function point_contain

    logical elemental function line_contain(this, that) result(result)
        class(line_t), intent(in)     :: this
        class(geometry_t), intent(in) :: that
        select type(that)
        type is (point_t)
            if (size(that%center) == 1) then
                result = (abs(this%center(1) - that%center(1)) <= this%length/2)
            else
                result = .false.
            end if
        type is (line_t)
            result = (abs(this%center(1) - that%center(1)) <= this%length/2 - that%length/2)
        class default
            result = .false.
        end select
    end function line_contain

    logical elemental function rectangle_contain(this, that) result(result)
        class(rectangle_t), intent(in) :: this
        class(geometry_t), intent(in)  :: that
        select type(that)
        type is (point_t)
            if (size(that%center) == 2) then
                result = (abs(this%center(1) - that%center(1)) <= this%length(1)/2) &
                    .and.(abs(this%center(2) - that%center(2)) <= this%length(2)/2)
            else
                result = .false.
            end if
        type is (line_t)
            result = (abs(this%center(1) - that%center(1)) <= this%length(1)/2 - that%length/2) &
                .and.(abs(this%center(2) - 0) <= this%length(2)/2)
        type is (rectangle_t)
            result = (abs(this%center(1) - that%center(1)) <= this%length(1)/2 - that%length(1)/2) &
                .and.(abs(this%center(2) - that%center(2)) <= this%length(2)/2 - that%length(2)/2)
        type is (circle_t)
            result = (abs(this%center(1) - that%center(1)) <= this%length(1)/2 - that%radius) &
                .and.(abs(this%center(2) - that%center(2)) <= this%length(2)/2 - that%radius)
        class default
            result = .false.
        end select
    end function rectangle_contain

    logical elemental function circle_contain(this, that) result(result)
        class(circle_t), intent(in)   :: this
        class(geometry_t), intent(in) :: that
        select type(that)
        type is (point_t)
            if (size(that%center) == 2) then
                result = (sqrt((this%center(1) - that%center(1))**2 &
                             + (this%center(2) - that%center(2))**2) <= this%radius)
            else
                result = .false.
            end if
        type is (line_t)
            result = (sqrt((this%center(1) - that%center(1))**2                  &
                         + (this%center(2) - 0)**2) <= this%radius) &
                .and.(sqrt((this%center(1) - that%center(1) - that%length)**2 &
                         + (this%center(2) - 0)**2) <= this%radius)
        type is (rectangle_t) !! Generate by Copilot, may be wrong
            result = (sqrt((this%center(1) - that%center(1))**2                  &
                         + (this%center(2) - that%center(2))**2) <= this%radius) &

                .and.(sqrt((this%center(1) - that%center(1) - that%length(1))**2 &
                         + (this%center(2) - that%center(2))**2) <= this%radius) &

                .and.(sqrt((this%center(1) - that%center(1) - that%length(1))**2                  &
                         + (this%center(2) - that%center(2) - that%length(2))**2) <= this%radius) &

                .and.(sqrt((this%center(1) - that%center(1))**2 &
                         + (this%center(2) - that%center(2) - that%length(2))**2) <= this%radius)
        type is (circle_t)
            result = (sqrt((this%center(1) - that%center(1))**2 &
                         + (this%center(2) - that%center(2))**2) <= this%radius - that%radius)
        class default
            result = .false.
        end select
    end function circle_contain

    logical elemental function cuboid_contain(this, that) result(result)
        class(cuboid_t), intent(in)   :: this
        class(geometry_t), intent(in) :: that
        select type(that)
        type is (point_t)
            if (size(that%center) == 3) then
                result = (abs(this%center(1) - that%center(1)) <= this%length(1)/2) &
                    .and.(abs(this%center(2) - that%center(2)) <= this%length(2)/2) &
                    .and.(abs(this%center(3) - that%center(3)) <= this%length(3)/2)
            else
                result = .false.
            end if
        type is (line_t)
            result = (abs(this%center(1) - that%center(1)) <= this%length(1)/2 - that%length/2) &
                .and.(abs(this%center(2) - 0) <= this%length(2)/2) &
                .and.(abs(this%center(3) - 0) <= this%length(3)/2)
        type is (rectangle_t)
            result = (abs(this%center(1) - that%center(1)) <= this%length(1)/2 - that%length(1)/2) &
                .and.(abs(this%center(2) - that%center(2)) <= this%length(2)/2 - that%length(2)/2) &
                .and.(abs(this%center(3) - 0) <= this%length(3)/2)
        type is (circle_t)
            result = (abs(this%center(1) - that%center(1)) <= this%length(1)/2 - that%radius) &
                .and.(abs(this%center(2) - that%center(2)) <= this%length(2)/2 - that%radius) &
                .and.(abs(this%center(3) - 0) <= this%length(3)/2)
        type is (cuboid_t)
            result = (abs(this%center(1) - that%center(1)) <= this%length(1)/2 - that%length(1)/2) &
                .and.(abs(this%center(2) - that%center(2)) <= this%length(2)/2 - that%length(2)/2) &
                .and.(abs(this%center(3) - that%center(3)) <= this%length(3)/2 - that%length(3)/2)
        type is (sphere_t)
            result = (abs(this%center(1) - that%center(1)) <= this%length(1)/2 - that%radius) &
                .and.(abs(this%center(2) - that%center(2)) <= this%length(2)/2 - that%radius) &
                .and.(abs(this%center(3) - that%center(3)) <= this%length(3)/2 - that%radius)
        class default
            result = .false.
        end select
    end function cuboid_contain

    logical elemental function sphere_contain(this, that) result(result)
        class(sphere_t), intent(in)   :: this
        class(geometry_t), intent(in) :: that
        select type(that)
        type is (point_t)
            if (size(that%center) == 3) then
                result = (sqrt((this%center(1) - that%center(1))**2 &
                             + (this%center(2) - that%center(2))**2 &
                             + (this%center(3) - that%center(3))**2) <= this%radius)
            else
                result = .false.
            end if
        type is (cuboid_t) !! Generate by Copilot, may be wrong
            result = (sqrt((this%center(1) - that%center(1))**2                  &
                         + (this%center(2) - that%center(2))**2                  &
                         + (this%center(3) - that%center(3))**2) <= this%radius) &

                .and.(sqrt((this%center(1) - that%center(1) - that%length(1))**2 &
                         + (this%center(2) - that%center(2))**2                  &
                         + (this%center(3) - that%center(3))**2) <= this%radius) &

                .and.(sqrt((this%center(1) - that%center(1) - that%length(1))**2 &
                         + (this%center(2) - that%center(2) - that%length(2))**2 &
                         + (this%center(3) - that%center(3))**2) <= this%radius) &

                .and.(sqrt((this%center(1) - that%center(1))**2                  &
                         + (this%center(2) - that%center(2) - that%length(2))**2 &
                         + (this%center(3) - that%center(3))**2) <= this%radius) &

                .and.(sqrt((this%center(1) - that%center(1))**2                                   &
                         + (this%center(2) - that%center(2))**2                                   &
                         + (this%center(3) - that%center(3) - that%length(3))**2) <= this%radius) &

                .and.(sqrt((this%center(1) - that%center(1) - that%length(1))**2                  &
                         + (this%center(2) - that%center(2))**2                                   &
                         + (this%center(3) - that%center(3) - that%length(3))**2) <= this%radius) &

                .and.(sqrt((this%center(1) - that%center(1) - that%length(1))**2                  &
                         + (this%center(2) - that%center(2) - that%length(2))**2                  &
                         + (this%center(3) - that%center(3) - that%length(3))**2) <= this%radius) &

                .and.(sqrt((this%center(1) - that%center(1))**2                  &
                         + (this%center(2) - that%center(2) - that%length(2))**2 &
                         + (this%center(3) - that%center(3) - that%length(3))**2) <= this%radius)

        type is (sphere_t)
            result = (sqrt((this%center(1) - that%center(1))**2 &
                         + (this%center(2) - that%center(2))**2 &
                         + (this%center(3) - that%center(3))**2) <= this%radius - that%radius)
        class default
            result = .false.
        end select
    end function sphere_contain

    logical elemental function point_intersect(this, that) result(result)
        class(point_t), intent(in)    :: this
        class(geometry_t), intent(in) :: that
        select type(that)
        type is (point_t)
            result = all(this%center == that%center)
        class default
            result = .false.
        end select
    end function point_intersect

    logical elemental function line_intersect(this, that) result(result)
        class(line_t), intent(in)     :: this
        class(geometry_t), intent(in) :: that
        select type(that)
        type is (line_t)
            result = abs(this%center(1) - that%center(1)) < this%length/2 + that%length/2
        class default
            result = .false.
        end select
    end function line_intersect

    logical elemental function rectangle_intersect(this, that) result(result)
        class(rectangle_t), intent(in) :: this
        class(geometry_t), intent(in)  :: that
        select type(that)
        type is (rectangle_t)
            result = (this%center(1) - this%length(1)/2) < (that%center(1) + that%length(1)/2) &
                 .or.(this%center(2) - this%length(2)/2) < (that%center(2) + that%length(2)/2) &
                 .or.(this%center(1) + this%length(1)/2) > (that%center(1) - that%length(1)/2) &
                 .or.(this%center(2) + this%length(2)/2) > (that%center(2) - that%length(2)/2)
        type is (circle_t)
            result = (abs(this%center(1) - that%center(1)) < this%length(1)/2 + that%radius) &
                .and.(abs(this%center(2) - that%center(2)) < this%length(2)/2 + that%radius)
        class default
            result = .false.
        end select
    end function rectangle_intersect

    logical elemental function circle_intersect(this, that) result(result)
        class(circle_t), intent(in)   :: this
        class(geometry_t), intent(in) :: that
        select type(that)
        type is (line_t)
            result = abs(this%center(1) - that%center(1)) < this%radius + that%length/2
        type is (rectangle_t)
            result = (abs(this%center(1) - that%center(1)) < this%radius + that%length(1)/2) &
                .and.(abs(this%center(2) - that%center(2)) < this%radius + that%length(2)/2)
        type is (circle_t)
            result = sqrt((this%center(1) - that%center(1))**2 &
                        + (this%center(2) - that%center(2))**2) < this%radius + that%radius
        class default
            result = .false.
        end select
    end function circle_intersect

    logical elemental function cuboid_intersect(this, that) result(result)
        class(cuboid_t), intent(in)   :: this
        class(geometry_t), intent(in) :: that
        select type(that)
        type is (line_t)
            result = (abs(this%center(1) - that%center(1)) < this%length(1)/2 + that%length/2) &
                .and.(abs(this%center(2) - 0) < this%length(2)/2 + that%length/2) &
                .and.(abs(this%center(3) - 0) < this%length(3)/2 + that%length/2)
        type is (rectangle_t)
            result = (abs(this%center(1) - that%center(1)) < this%length(1)/2 + that%length(1)/2) &
                .and.(abs(this%center(2) - that%center(2)) < this%length(2)/2 + that%length(2)/2) &
                .and.(abs(this%center(3)) < this%length(3)/2)
        type is (circle_t)
            result = (abs(this%center(1) - that%center(1)) < this%length(1)/2 + that%radius) &
                .and.(abs(this%center(2) - that%center(2)) < this%length(2)/2 + that%radius) &
                .and.(abs(this%center(3)) < this%length(3)/2)
        type is (sphere_t)
            result = (abs(this%center(1) - that%center(1)) < this%length(1)/2 + that%radius) &
                .and.(abs(this%center(2) - that%center(2)) < this%length(2)/2 + that%radius) &
                .and.(abs(this%center(3) - that%center(3)) < this%length(3)/2 + that%radius)
        class default
            result = .false.
        end select
    end function cuboid_intersect

    logical elemental function sphere_intersect(this, that) result(result)
        class(sphere_t), intent(in)   :: this
        class(geometry_t), intent(in) :: that
        select type(that)
        type is (line_t)
            result = abs(this%center(1) - that%center(1)) < this%radius + that%length/2
        type is (rectangle_t)
            result = (abs(this%center(1) - that%center(1)) < this%radius + that%length(1)/2) &
                .and.(abs(this%center(2) - that%center(2)) < this%radius + that%length(2)/2) &
                .and.(abs(this%center(3)) < this%radius)
        type is (circle_t)
            result = sqrt((this%center(1) - that%center(1))**2 &
                        + (this%center(2) - that%center(2))**2) < this%radius + that%radius
        type is (sphere_t)
            result = sqrt((this%center(1) - that%center(1))**2 &
                        + (this%center(2) - that%center(2))**2 &
                        + (this%center(3) - that%center(3))**2) < this%radius + that%radius
        class default
            result = .false.
        end select
    end function sphere_intersect

    ! pure subroutine point_clean(this)
    !     class(point_t), intent(inout) :: this
    !     deallocate(this%center)
    ! end subroutine point_clean

    ! pure subroutine line_clean(this)
    !     class(line_t), intent(inout) :: this
    !     ! deallocate(this%center)
    ! end subroutine line_clean

    ! pure subroutine rectangle_clean(this)
    !     class(rectangle_t), intent(inout) :: this
    !     ! deallocate(this%center)
    ! end subroutine rectangle_clean

    ! pure subroutine circle_clean(this)
    !     class(circle_t), intent(inout) :: this
    !     ! deallocate(this%center)
    ! end subroutine circle_clean

    ! pure subroutine cuboid_clean(this)
    !     class(cuboid_t), intent(inout) :: this
    !     ! deallocate(this%center)
    ! end subroutine cuboid_clean

    ! pure subroutine sphere_clean(this)
    !     class(sphere_t), intent(inout) :: this
    !     ! deallocate(this%center)
    ! end subroutine sphere_clean

end module geometry_m