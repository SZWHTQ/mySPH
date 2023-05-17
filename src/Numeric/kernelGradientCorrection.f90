module KGC_m
    use SPH
    use ctrl_dict, only: Field
    use tools_m, only: MCPE
    implicit none

contains
    subroutine kernelGradientCorrection(ntotal, P)
        integer, intent(in) :: ntotal
        type(Particle), intent(inout) :: P(:)
        real(8) :: dx(Field%Dim)
        real(8) :: Matrix(Field%Dim, Field%Dim), Temp(Field%Dim, Field%Dim)
        real(8) :: correctedKernelGradient(Field%Dim)

        integer i, j, k, d, e

        do i = 1, ntotal
            Matrix = 0
            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)
                dx(:) = P(j)%x(:) - P(i)%X(:)
                do d = 1, Field%Dim
                    do e = 1, Field%Dim
                        Matrix(d, e) = Matrix(d, e) &
                            + dx(e) * P(i)%dwdx(d, k) * P(j)%Mass/P(j)%Density
                    end do
                end do
            end do
            do k = 1, P(i)%neighborNum
                j = P(i)%neighborList(k)
                Temp = Matrix
                call MCPE(Matrix, P(i)%dwdx(:, k), correctedKernelGradient)
                Matrix = Temp
                P(i)%dwdx(:, k) = correctedKernelGradient(:)
            end do
        end do

    end subroutine kernelGradientCorrection
    
end module KGC_m