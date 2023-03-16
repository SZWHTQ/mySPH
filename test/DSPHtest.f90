program main
    use ctrl_dict
    use nnps_m
    use kernel_m
    implicit none
    type :: sph_t
        real(8) :: y(0:1)
    end type sph_t
    type :: csph_t
        real(8) :: y(0:1)
        real(8) :: sum(0:1)
    end type csph_t
    type :: dsph_t
        real(8) :: y(0:1)
        real(8) :: sum(0:1)
        integer :: d !! Discontinuous location
    end type dsph_t
    interface
        elemental function func(var) result(retval)
            implicit none
            real(8), intent(in) :: var
            real(8) :: retval
        end function func
        elemental function diff_func(var) result(retval)
            implicit none
            real(8), intent(in) :: var
            real(8) :: retval
        end function diff_func
    end interface
    type(sph_t),  allocatable :: sph(:)
    type(csph_t), allocatable :: csph(:)
    type(dsph_t), allocatable :: dsph(:)
    real(8), allocatable :: exact(:,:)
    real(8) :: self, hv(1)
    real(8) :: f_max, f_min, criteria, ratio !! Discontinuous criteria
    real(8) :: delta
    integer, allocatable :: neighborNum(:)
    integer i, j, k, s, step
    

    dim  = 1
    nnps = 1
    skf  = 1
    maxn = 50
    ntotal = maxn
    max_interaction = maxn * 5

    call initialize()
    allocate(exact(0:1, ntotal))
    allocate(  sph(ntotal) )
    allocate( csph(ntotal) )
    allocate( dsph(ntotal) )
    allocate(neighborNum(ntotal))

    delta = 1. / ntotal
    hsml = delta * 2
    do i = 1, ntotal
        x(1, i) = -1. + dble(2*i-1) / ntotal
    end do

    do i = 1, ntotal
        exact(0,i) = func(x(1, i))
        exact(1,i) = diff_func(x(1, i))
    end do
    
    ! write(*,"(5(E15.3))") exact

    niac = 0
    w    = 0
    dwdx = 0
    call search_particles(nnps, ntotal+ndummy, x, hsml, max_interaction, &
    niac, neighborList, w, dwdx, neighborNum)

    ! do i = 1, max_interaction
    !     write(*, "(2(I0, 2X), 2(G0, 2X))") neighborList(i,:), w(i), dwdx(:, i)
    ! end do

    criteria = 0.19
    f_max  = maxval(exact(0,:))
    f_min  = minval(exact(0,:))
    dsph%d = ntotal + 1
    s = 0
    do k = 1, niac
        i = neighborList(k, 1)
        j = neighborList(k, 2)
        if ( abs((exact(0,j) - exact(0,i))/(f_max-f_min)) >= criteria &
       .and. abs((exact(0,j) - exact(0,i))/(f_max-f_min)) >= ratio ) then
            ratio     = abs((exact(0,j) - exact(0,i))/(f_max-f_min))
            if ( s /= i ) then
                dsph(i)%d = j
                s = i
            end if
        end if
    end do
    ! write(*,"(8(I0, 2X))") dsph%d

    do i = 1, ntotal
        call kernel(dble(0), [dble(0)], hsml(i), self, hv(:))
        sph(i)%y(0) = exact(0,i) * self * delta
        sph(i)%y(1) = 0
        
        csph(i)%y(0)   = exact(0,i) * self * delta
        csph(i)%y(1)   = 0
        csph(i)%sum(0) = self * delta
        csph(i)%sum(1) = 0

        dsph(i)%y(0)   = exact(0,i) * self * delta
        dsph(i)%y(1)   = 0
        dsph(i)%sum(0) = self * delta
        dsph(i)%sum(1) = 0
    end do

    do k = 1, niac
        i = neighborList(k, 1)
        j = neighborList(k, 2)
        sph(i)%y(0) = sph(i)%y(0) + exact(0,j) * w(k) * delta
        sph(j)%y(0) = sph(j)%y(0) + exact(0,i) * w(k) * delta
        sph(i)%y(1) = sph(i)%y(1) + (exact(0,j)-exact(0,i)) * dwdx(1, k) * delta
        sph(j)%y(1) = sph(j)%y(1) + (exact(0,j)-exact(0,i)) * dwdx(1, k) * delta

        csph(i)%y(0)   = csph(i)%y(0) + exact(0,j) * w(k) * delta
        csph(j)%y(0)   = csph(j)%y(0) + exact(0,i) * w(k) * delta
        csph(i)%y(1)   = csph(i)%y(1) + (exact(0,j)-exact(0,i)) * dwdx(1, k) * delta
        csph(j)%y(1)   = csph(j)%y(1) + (exact(0,j)-exact(0,i)) * dwdx(1, k) * delta
        csph(i)%sum(0) = csph(i)%sum(0) + w(k) * delta
        csph(j)%sum(0) = csph(j)%sum(0) + w(k) * delta
        csph(i)%sum(1) = csph(i)%sum(1) + (x(1,j)-x(1,i)) * dwdx(1, k) * delta
        csph(j)%sum(1) = csph(j)%sum(1) + (x(1,j)-x(1,i)) * dwdx(1, k) * delta

        dsph(i)%y(0) = dsph(i)%y(0) + exact(0,j) * w(k) * delta
        dsph(j)%y(0) = dsph(j)%y(0) + exact(0,i) * w(k) * delta
        dsph(i)%y(1) = dsph(i)%y(1) + (exact(0,j)-exact(0,i)) * dwdx(1, k) * delta
        dsph(j)%y(1) = dsph(j)%y(1) + (exact(0,j)-exact(0,i)) * dwdx(1, k) * delta
        if ( j >= dsph(i)%d ) then
            dsph(i)%y(0) = dsph(i)%y(0) &
                         - (exact(0,dsph(i)%d)-exact(0,i))* w(k) * delta
            dsph(i)%y(1) = dsph(i)%y(1) &
                         - (exact(0,dsph(i)%d)-exact(0,i))* dwdx(1, k) * delta &
                         - ((x(1,j)-x(1,dsph(i)%d))*exact(1,dsph(i)%d) - (x(1,j)-x(1,i))*exact(1,i)) * dwdx(1, k) * delta
            dsph(j)%y(0) = dsph(j)%y(0) &
                         - (exact(0,i)-exact(0,j))* w(k) * delta
            dsph(j)%y(1) = dsph(j)%y(1) &
                         - (exact(0,j)-exact(0,i))* dwdx(1, k) * delta &
                         - ((x(1,i)-x(1,i))*exact(1,i) - (x(1,j)-x(1,i))*exact(1,j)) * dwdx(1, k) * delta
        end if
        dsph(i)%sum(0) = dsph(i)%sum(0) + w(k) * delta
        dsph(j)%sum(0) = dsph(j)%sum(0) + w(k) * delta
        dsph(i)%sum(1) = dsph(i)%sum(1) + (x(1,j)-x(1,i)) * dwdx(1, k) * delta
        dsph(j)%sum(1) = dsph(j)%sum(1) + (x(1,j)-x(1,i)) * dwdx(1, k) * delta
    end do

    do i = 1, ntotal
        csph(i)%y(0) = csph(i)%y(0) / csph(i)%sum(0)
        csph(i)%y(1) = csph(i)%y(1) / csph(i)%sum(1)
        dsph(i)%y(0) = dsph(i)%y(0) / dsph(i)%sum(0)
        dsph(i)%y(1) = dsph(i)%y(1) / dsph(i)%sum(1)
    end do
    ! write(*, "(8(E15.7, 2X))") csph%sum(0)
    ! write(*,*)
    ! write(*, "(8(E15.7, 2X))") csph%sum(1)

    step = 0
    do s = 1, step
        do i = 1, ntotal
            call kernel(dble(0), [dble(0)], hsml(i), self, hv(:))
            sph(i)%y(0) = sph(i)%y(0) * self * delta
            sph(i)%y(1) = 0
                  
            csph(i)%y(0)   = csph(i)%y(0) * self * delta
            csph(i)%y(1)   = 0
            csph(i)%sum(0) = self * delta
            csph(i)%sum(1) = 0

            dsph(i)%y(0)   = dsph(i)%y(0) * self * delta
            dsph(i)%y(1)   = 0
            dsph(i)%sum(0) = self * delta
            dsph(i)%sum(1) = 0
        end do
        do k = 1, niac
            i = neighborList(k, 1)
            j = neighborList(k, 2)
            sph(i)%y(0) = sph(i)%y(0) + sph(j)%y(0) * w(k) * delta
            sph(j)%y(0) = sph(j)%y(0) + sph(i)%y(0) * w(k) * delta
            sph(i)%y(1) = sph(i)%y(1) + (sph(j)%y(0)-sph(i)%y(0)) * dwdx(1, k) * delta
            sph(j)%y(1) = sph(j)%y(1) + (sph(j)%y(0)-sph(i)%y(0)) * dwdx(1, k) * delta
            
            csph(i)%y(0)   = csph(i)%y(0) + csph(j)%y(0) * w(k) * delta
            csph(j)%y(0)   = csph(j)%y(0) + csph(j)%y(0) * w(k) * delta
            csph(i)%y(1)   = csph(i)%y(1) + (csph(j)%y(0) - csph(i)%y(0)) * dwdx(1, k) * delta
            csph(j)%y(1)   = csph(j)%y(1) + (csph(j)%y(0) - csph(i)%y(0)) * dwdx(1, k) * delta
            csph(i)%sum(0) = csph(i)%sum(0) + w(k) * delta
            csph(j)%sum(0) = csph(j)%sum(0) + w(k) * delta
            csph(i)%sum(1) = csph(i)%sum(1) + (x(1, j) - x(1, i)) * dwdx(1, k) * delta
            csph(j)%sum(1) = csph(j)%sum(1) + (x(1, j) - x(1, i)) * dwdx(1, k) * delta
        end do
        do i = 1, ntotal
            csph(i)%y(0) = csph(i)%y(0) / csph(i)%sum(0)
            csph(i)%y(1) = csph(i)%y(1) / csph(i)%sum(1)
        end do
    end do


    open(unit=11, file="SPH.dat")
    write(11, "(A7, 7(A15))") "Index", "X", "Exact", "Y", "Error", "Diff_Exact", "Diff_Y", "Diff_Error"
    do i = 1, ntotal
        write(11, "(I8, 7(E15.7))") i, x(1, i), &
                                    exact(0,i), sph(i)%y(0), abs((exact(0,i)-sph(i)%y(0))/exact(0,i)), &
                                    exact(1,i), sph(i)%y(1), abs((exact(1,i)-sph(i)%y(1))/exact(1,i))
    end do

    open(unit=12, file="CSPH.dat")
    write(12, "(A7, 7(A15))") "Index", "X", "Exact", "Y", "Error", "Diff_Exact", "Diff_Y", "Diff_Error"
    do i = 1, ntotal
        write(12, "(I8, 7(E15.7))") i, x(1, i), &
                                    exact(0,i), csph(i)%y(0), abs((exact(0,i)-csph(i)%y(0))/exact(0,i)), &
                                    exact(1,i), csph(i)%y(1), abs((exact(1,i)-csph(i)%y(1))/exact(1,i))
    end do

    open(unit=13, file="DSPH.dat")
    write(13, "(A7, 7(A15))") "Index", "X", "Exact", "Y", "Error", "Diff_Exact", "Diff_Y", "Diff_Error"
    do i = 1, ntotal
        write(13, "(I8, 7(E15.7))") i, x(1, i), &
                                    exact(0,i), dsph(i)%y(0), abs((exact(0,i)-dsph(i)%y(0))/exact(0,i)), &
                                    exact(1,i), dsph(i)%y(1), abs((exact(1,i)-dsph(i)%y(1))/exact(1,i))
    end do

    deallocate(exact, sph, csph, dsph)
    
end program main

elemental function func(var) result(retval)
    implicit none
    real(8), parameter :: PI = acos(-1.)
    real(8), intent(in) :: var
    real(8) :: retval

    if ( -1 <= var .and. var <= 0.5) then
        retval = exp(3 * var) - 1
    else if ( 0.5 < var .and. var <= 1) then
        retval = 1 - log(2 * var)
    else
        retval = 0
    end if

    ! retval = sin(2*PI*var)
    
end function func

elemental function diff_func(var) result(retval)
    implicit none
    real(8), parameter :: PI = acos(-1.)
    real(8), intent(in) :: var
    real(8) :: retval

    if ( -1 <= var .and. var <= 0.5) then
        retval = 3 * exp(3 * var)
    else if ( 0.5 < var .and. var <= 1) then
        retval = - 1 / var
    else
        retval = 0
    end if

    ! retval = cos(2*PI*var)
    
end function diff_func