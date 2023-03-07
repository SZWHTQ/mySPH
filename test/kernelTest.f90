program main
    use ctrl_dict, only: skf, dim
    use kernel_m
    use tools_m, only: PI
    implicit none
    type :: sph_t
        real(8) :: y(0:1)
    end type sph_t
    type :: cpsh_t
        real(8) :: y(0:1)
        real(8) :: sum(0:1)
    end type cpsh_t
    type :: dpsh_t
        real(8) :: y(0:1)
        real(8) :: sum(0:1)
        integer :: d !! Discontinuous location
    end type dpsh_t
    integer, parameter :: ntotal = 21
    real(8) :: x(ntotal), y(ntotal)
    real(8) :: w(ntotal), dwdx(1,ntotal)
    real(8) :: delta, dx(1), hsml
    type(sph_t)  :: sph
    type(cpsh_t) :: csph
    integer i, k, m

    
    write(*,"(A5, 10X, 10(A17))") "Index", "Exact", "Exact diff", &
    "SPH",  "SPH diff",  "SPH Error",  "SPH diff Error", &
   "CSPH", "CSPH diff", "CSPH Error", "CSPH diff Error"
    do k = 1, 30
        dim = size(dx)
        delta = 1./k
        hsml  = delta * 1
        m = (ntotal+1)/2
        do i = 1, ntotal
            x(i) = delta * i + 0.5
            y(i) = func(x(i))
        end do
            
        skf = 1
        do i = 1, ntotal
            dx = x(m) - x(i)
            call kernel(norm2(dx), dx, hsml, w(i), dwdx(:, i))
        end do
    
        sph%y(:) = 0
        csph%y(:)   = 0
        csph%sum(:) = 0
        do i = 1, ntotal
            sph%y(0) = sph%y(0) + y(i) * w(i) * delta
            sph%y(1) = sph%y(1) + y(i) * dwdx(1,i) * delta

            csph%y(0)   = csph%y(0) + y(i) * w(i) * delta
            csph%y(1)   = csph%y(1) + y(i) * dwdx(1,i) * delta
            csph%sum(0) = csph%sum(0) + w(i) * delta
            csph%sum(1) = csph%sum(1) + (x(i)-x(m))*dwdx(1,i) * delta
        end do
        ! do i = 1, ntotal
        !     csph%y(0) = csph%y(0) / csph%sum(0)
        !     csph%y(1) = csph%y(1) / csph%sum(1)
        ! end do
        write(*, "(I3, 12X, 2(E15.7, 2X))",     advance="no") k, func(x(m)), diff_func(x(m))
        write(*, "(2(E15.7, 2X), 2(F15.4, A, 2X))", advance="no")  sph%y(0),  sph%y(1), &
            100*abs((func(x(m)) -  sph%y(0)) / func(x(m))), "%", &
            100*abs((diff_func(x(m)) -  sph%y(1)) / diff_func(x(m))), "%"
        write(*, "(2(E15.7, 2X), 2(F15.4, A, 2X))") csph%y(0), csph%y(1), &
            100*abs((func(x(m)) - csph%y(0)) / func(x(m))), "%", &
            100*abs((diff_func(x(m)) - csph%y(1)) / diff_func(x(m))), "%"

        ! write(*, "(2(G0, 2X), F0.4, A)") sph%y(1) , diff_func(x(m)), &
        !                        100*abs(diff_func(x(m)) - sph%y(1) ) / diff_func(x(m)), "%"
        ! do i = 1, ntotal
        !     write(*, "(I2, 2X, G0, 2X, G0)") i, w(i), dwdx(:,i)
        ! end do
    end do
    

contains
    function func(var) result(re)
        real(8) var, re
        re = cos(var)
    end function
    function diff_func(var) result(re)
        real(8) var, re
        re = -sin(var)
    end function
end program main