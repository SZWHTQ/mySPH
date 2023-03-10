module density_m
    use ctrl_dict, only: dim, norm_dens_w, DSPH_w
    ! use initial_m, only: v, mass, rho, hsml, neighborList, w, dwdx
    use kernel_m,  only: kernel

    implicit none

contains
    !!! Subroutine to calculate the density with SPH summation algorithm
    subroutine sum_density(ntotal, mass, rho, hsml, neighborNum, neighborList, w, div_r)
        integer, intent(in)    :: ntotal
        real(8), intent(in)    :: mass(:)        !! Mass of each particle
        real(8), intent(in)    :: hsml(:)        !! Smoothing length of each particle
        integer, intent(in)    :: neighborNum(:) !! Number of neighbors of each particle
        integer, intent(in)    :: neighborList(:,:)      !! Pair of particles
        real(8), intent(in)    :: w(:, :)        !! Kernel value of each neighborList
        real(8), intent(in)    :: div_r(:)     !! Divergence of each particle
        real(8), intent(inout) :: rho(:)       !! Density of each particle
        real(8) :: self
        real(8) :: hv(dim)
        real(8), allocatable :: wi(:)  !! Integration of the kernel itself
        integer i, j, k

        allocate(wi(ntotal), source=0._8)
        hv = 0

        !!! Self density of each particle: Wii (Kernel for distance 0)
        !!! and take contribution of particle itself

        !!! Firstly, calculate the integration of the kernel over the space
        if ( norm_dens_w ) then
            !$OMP PARALLEL DO PRIVATE(i, self, hv)
            do i = 1, ntotal
                if ( div_r(i) < 1.5 ) then
                    call kernel(dble(0), 1*hv, hsml(i), self, hv)
                    wi(i) = mass(i)/rho(i) * self
                end if
            end do
            !$OMP END PARALLEL DO
        end if

        if ( norm_dens_w ) then
            !$OMP PARALLEL DO PRIVATE(i, j, k) REDUCTION(+:wi)
            do i = 1, ntotal
                if ( div_r(i) < 1.5 ) then
                    do k = 1, neighborNum(i)
                        j = neighborList(i, k)
                        wi(i) = wi(i) + mass(j)/rho(j) * w(i, k)
                    end do
                end if
            end do
            !$OMP END PARALLEL DO
        end if

        !!! Secondly, calculate the rho integration over the space
        !$OMP PARALLEL DO PRIVATE(i, self, hv)
        do i = 1, ntotal
            call kernel(dble(0), 1*hv, hsml(i), self, hv)
            rho(i) = mass(i) * self
        end do
        !$OMP END PARALLEL DO

        !!! Calculate SPH sum for rho:
        ! !$OMP PARALLEL DO PRIVATE(i, j, k) REDUCTION(+:rho)
        do i = 1, ntotal
            do k = 1, neighborNum(i)
                j = neighborList(i, k)
                rho(i) = rho(i) + mass(j) * w(i, k)
            end do
        end do
        ! !$OMP END PARALLEL DO

        !!! Thirdly, calculate the normalized rho, rho = Σrho / Σw
        if ( norm_dens_w ) then
            !$OMP PARALLEL DO PRIVATE(i)
            do i = 1, ntotal
                if ( div_r(i) < 1.5 ) then
                    rho(i) = rho(i) / wi(i)
                end if
            end do
            !$OMP END PARALLEL DO
        end if

        deallocate(wi)

    end subroutine sum_density

    !!! Subroutine to calculate the density with SPH continuity approach
    subroutine con_density(ntotal, v, mass, neighborNum, neighborList, dwdx, drhodt)
        integer, intent(in)  :: ntotal
        real(8), intent(in)  :: v(:, :)        !! Velocity of each particle
        real(8), intent(in)  :: mass(:)        !! Mass of each particle
        integer, intent(in)  :: neighborNum(:) !! Number of neighbors of each particle
        integer, intent(in)  :: neighborList(:,:)      !! Pair of particles
        real(8), intent(in)  :: dwdx(:, :, :)  !! Derivative of kernel value of each neighborList
        real(8), intent(inout) :: drhodt(:)    !! Density change rate of each particle

        integer i, j, k

        do i = 1, ntotal
            drhodt(i) = 0
        end do

        !$OMP PARALLEL DO PRIVATE(i, j, k) REDUCTION(+:drhodt)
        do i = 1, ntotal
            do k = 1, neighborNum(i)
                j = neighborList(i, k)
                drhodt(i) = drhodt(i) &
                    + mass(j)         &
                    * sum((v(:, i) - v(:, j)) * dwdx(:, i, k))
            end do
        end do
        !$OMP END PARALLEL DO


    end subroutine con_density

    !!! PVRS Riemann Solver
    subroutine con_density_riemann(ntotal, x, v, mass, rho, p, c, neighborNum, neighborList, dwdx, drhodt)
        integer, intent(in)  :: ntotal
        real(8), intent(in)  :: x(:, :)        !! Position of each particle
        real(8), intent(in)  :: v(:, :)        !! Velocity of each particle
        real(8), intent(in)  :: mass(:)        !! Mass of each particle
        real(8), intent(in)  :: rho(:)         !! Density of each particle
        real(8), intent(in)  :: p(:)           !! Pressure of each particle
        real(8), intent(in)  :: c(:)           !! Sound speed of each particle
        integer, intent(in)  :: neighborNum(:) !! Number of neighbors of each particle
        integer, intent(in)  :: neighborList(:,:)      !! Pair of particles
        real(8), intent(in)  :: dwdx(:, :, :)  !! Derivative of kernel value of each neighborList
        real(8), intent(inout) :: drhodt(:)    !! Density change rate of each particle
        real(8) :: Z_l, Z_r, v_l, v_r
        real(8) :: v_ij
        real(8) :: v_star(dim), e_ij(dim)
        integer i, j, k

        do i = 1, ntotal
            drhodt(i) = 0
        end do
        v_ij = 0
        e_ij = 0

        do i = 1, ntotal
            do k = 1, neighborNum(i)
                j = neighborList(i, k)

                Z_l = rho(i) * c(i)
                Z_r = rho(j) * c(j)

                e_ij = (x(:, j) - x(:, i)) / norm2(x(:, j) - x(:, i))

                v_l = dot_product(v(:, i), e_ij)
                v_r = dot_product(v(:, j), e_ij)

                v_ij = ( Z_l*v_l + Z_r*v_r + (p(i)-p(j)) ) &
                     / ( Z_l + Z_r )
                v_star = v_ij*e_ij + ((v(:, i)+v(:, j))/2 - ((v_l+v_r)/2)*e_ij)

                drhodt(i) = drhodt(i) &
                + 2*rho(i) * mass(j)/rho(j) * dot_product((v(:, i) - v_star), dwdx(:, i, k))
            end do
        end do

    end subroutine con_density_riemann

    subroutine sum_density_dsph(ntotal, mass, rho, hsml, neighborNum, neighborList, w)
        integer, intent(in)  :: ntotal
        real(8), intent(in)  :: mass(:)        !! Mass of each particle
        real(8), intent(in)  :: hsml(:)        !! Smoothing length of each particle
        integer, intent(in)  :: neighborNum(:) !! Number of neighbors of each particle
        integer, intent(in)  :: neighborList(:,:)      !! Pair of particles
        real(8), intent(in)  :: w(:, :)        !! Kernel value of each neighborList
        real(8), intent(inout) :: rho(:)       !! Density of each particle
        real(8) :: self
        real(8) :: hv(dim)
        real(8) :: rho_max, rho_min, criteria, ratio
        real(8), allocatable :: wi(:)  !! Integration of the kernel itself
        integer, allocatable :: dc_point(:)
        integer i, j, k, s

        allocate(wi(ntotal), source=0._8)
        allocate(dc_point(ntotal), source=0)
        hv = 0

        !!! Self density of each particle: Wii (Kernel for distance 0)
        !!! and take contribution of particle itself

        !!! Firstly, calculate the integration of the kernel over the space
        do i = 1, ntotal
            call kernel(dble(0), 1*hv, hsml(i), self, hv)
            wi(i) = mass(i)/rho(i) * self
        end do

        criteria = 0.2
        ratio    = 0
        rho_max  = maxval(rho)
        rho_min  = minval(rho)
        dc_point = ntotal + 1
        s = 0

        do i = 1, ntotal
            do k = 1, neighborNum(i)
                j = neighborList(i, k)
                wi(i) = wi(i) + mass(j)/rho(j) * w(i, k)
                if ( abs((rho(j)-rho(i))/(rho_max-rho_min)) >= criteria &
               .and. abs((rho(j)-rho(i))/(rho_max-rho_min)) >= ratio ) then
                    ratio = abs((rho(j)-rho(i))/(rho_max-rho_min))
                    if ( s /= i ) then
                          dc_point(i) = j
                          s = i
                    end if
                 end if
            end do
        end do

        !!! Secondly, calculate the rho integration over the space
        do i = 1, ntotal
            call kernel(dble(0), 1*hv, hsml(i), self, hv)
            rho(i) = mass(i) * self
        end do

        do i = 1, ntotal
            do k = 1, neighborNum(i)
                j = neighborList(i, k)
                s = dc_point(i)
                rho(i) = rho(i) + mass(j)*w(i, k)
                if ( j >= s ) then
                    rho(i) = rho(i) - (rho(s)-rho(i)) * mass(j)/rho(j) * w(i, k)
                end if
            end do
        end do


        !!! Thirdly, calculate the normalized rho, rho = Σrho / Σw
        do i = 1, ntotal
            rho(i) = rho(i) / wi(i)
        end do

        deallocate(wi, dc_point)

    end subroutine sum_density_dsph

end module density_m