module KGC_m
    use SPH
    use ctrl_dict, only: Field
    use tools_m, only: MCPE
    implicit none

contains
    subroutine kernelGradientCorrection(ntotal, P)
        integer, intent(in) :: ntotal
        type(Particle), intent(inout) :: P(:)
        real(8), allocatable :: dx(:)
        real(8), allocatable :: A(:, :, :), Temp(:, :), B(:)
        integer, allocatable :: v(:)
        integer :: iflag

        integer i, j, k, d, e

        allocate(dx(Field%Dim), source=0._8)
        allocate(A(Field%Dim, Field%Dim, ntotal), Temp(Field%Dim, Field%Dim), source=0._8)
        allocate(B(Field%Dim), source=0._8)
        allocate(v(Field%Dim), source=0)

        !$omp parallel do private(i, j, k, d, e, dx, Temp, B) reduction(+:A)
        do i = 1, ntotal
            ! if ( P(i)%Type <= 100 ) cycle
            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)
                ! if ( (abs(P(i)%Type)-100)*(abs(P(j)%Type)-100) < 0 ) cycle
                ! if ( P(i)%Type /= P(j)%Type ) then
                !     cycle
                ! end if
                if ( P(i)%Type /= P(j)%Type ) cycle
                dx(:) = P(j)%x(:) - P(i)%X(:)
                do d = 1, Field%Dim
                    do e = 1, Field%Dim
                        A(d, e, i) = A(d, e, i) &
                            + dx(e) * P(i)%dwdx(d, k) * P(j)%Mass/P(j)%Density
                    end do
                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, k, Temp, B, v, iflag)
        do i = 1, ntotal
            do k = 1, P(i)%neighborNum
                v = 0
                iflag = 0
                Temp = A(:, :, i)
                B = P(i)%dwdx(:, k)
                !! MCPE
                call MCPE(Temp, B, P(i)%dwdx(:, k)) !! Temp will change with call MCPE
                ! !! Lapack
                ! call DGESV(Field%Dim, 1, Temp, &
                !            Field%Dim, v, B,    &
                !            Field%Dim, iflag)
                ! if ( iflag < 0 ) then
                !     write(*, "(I0)") iflag
                !     error stop "the i-th argument had an illegal value"
                ! else if ( iflag > 0 ) then
                !     write(*, "(I0)") iflag
                !     error stop "U(i,i) is exactly zero"
                ! end if
                ! P(i)%dwdx(:, k) = B(:)
            end do
        end do
        !$omp end parallel do

        deallocate(dx)
        deallocate(A, Temp)
        deallocate(B)

    end subroutine kernelGradientCorrection
    
end module KGC_m