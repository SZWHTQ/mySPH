program main
    use eos_m, only: mie_gruneisen_eos_of_solid, mie_gruneisen_eos_of_water, water_polynomial_eos
    implicit none
    real(8) :: rho, e, p

    call mie_gruneisen_eos_of_solid(3d4, 0._8, p)
    write(*,"(A, F0.3, 2A, F0.3, A)") ">> Pressure: ", p*1e-6, " MPa"

    write(*,"(A)", advance="no") "Input density (kg/m^3) and internal energy (J/kg): "
    read(*,*) rho, e

    call mie_gruneisen_eos_of_water(rho, e, p)
    write(*,"(A, F0.3, 2A, F0.3, A)") ">> Pressure: ", p*1e-6, " MPa", &
                                        "   Depth: ", p / ((rho+1000)/2) / 9.81, " m"

    call water_polynomial_eos(rho, e, p)
    write(*,"(A, F0.3, 2A, F0.3, A)") ">> Pressure: ", p*1e-6, " MPa", &
                                        "   Depth: ", p / ((rho+1000)/2) / 9.81, " m"

    
end program main