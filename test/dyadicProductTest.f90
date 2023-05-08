program dyadicProductTest
    use tools_m, only: dyadic_product
    implicit none
    real(8) :: a(4) = [1.0, 2.0, 3.0, 4.0]
    real(8) :: b(3) = [5.0, 6.0, 7.0]
    real(8) :: c(4, 3)

    integer i, j

    c = dyadic_product(a, b)

    write(*, "(3(G0, 2X))") ((c(i,j), j=1,3), i=1,4)

    
end program dyadicProductTest