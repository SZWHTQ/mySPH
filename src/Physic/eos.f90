module eos_m
    ! use tools_m

    implicit none

contains
    !!! Gamma law EOS: subroutine to calculate the pressure and sound
    !!! Type = 1
    elemental subroutine gas_eos(rho, e, p, c)
        real(8), intent(in)  :: rho, e
        real(8), intent(inout) :: p,   c
        real(8), parameter :: gamma = 1.4

        p = (gamma-1) * rho * e
        ! c = sqrt((gamma-1)*e)
        c = sqrt(gamma*p/rho)

    end subroutine gas_eos

    !!! Artificial equation of state for the artificial compressibility
    !!! Artificial EOS, Form 1 (Monaghan, 1994) !! Modified for dam break
    !!! Type = 2
    elemental subroutine arti_water_eos_1(rho, p, Density)
        use ctrl_dict, only: Project
        real(8), intent(inout) :: rho
        real(8), intent(inout) :: p
        logical, intent(in), optional :: Density
        real(8) :: c !! Artificial/Lagrangian sound speed
        real(8), parameter :: gamma = 7
        real(8), parameter :: rho0 = 1000
        real(8) :: b

        ! c = 10
        ! c = 50 !!for "Two-dimensional dam break"
        select case(Project%nick)
        case ("dam_break")
            c = 50
        case ("water_impact")
            c = 10
        end select
        b = c**2 * rho0 / gamma

        if ( .not. present(Density) ) then
            p = b * ((rho/rho0)**gamma - 1)
        else
            if ( Density ) then
                rho = (p / b + 1)**(1._8/gamma) * rho0
            else
                p = b * ((rho/rho0)**gamma - 1)
            end if
        end if

    end subroutine arti_water_eos_1

    !!! Artificial EOS, form 2 (Morris, 1997)
    !!! Type = 3
    elemental subroutine arti_water_eos_2(rho, p, c)
        real(8), intent(in)  :: rho
        real(8), intent(inout) :: p, c

        c = 0.01
        p = c**2 * rho

    end subroutine arti_water_eos_2

    !!! TNT gas EOS/Gamma Law
    !!! Type = 4
    elemental subroutine tnt_eos(rho, e, p)
        real(8), intent(in)  :: rho, e
        real(8), intent(inout) :: p
        real(8), parameter :: gamma = 1.4

        p = (gamma-1) * rho * e

    end subroutine tnt_eos

    !!! TNT gas EOS/Jones-Wilkins-Lee
    !!! Type = 5
    elemental subroutine jwl_eos(rho, e, p)
        real(8), intent(in)  :: rho, e
        real(8), intent(inout) :: p
        real(8) :: ratio
        real(8), parameter :: rho0 = 1630
        real(8), parameter :: A  = 3.712e11, B  = 3.21e9
        real(8), parameter :: R1 = 4.15,     R2 = 0.95
        real(8), parameter :: omega = 0.3

        ratio = rho/rho0

        p = A*(1-omega*ratio/R1)*exp(-R1/ratio) &
          + B*(1-omega*ratio/R2)*exp(-R2/ratio) &
          + omega * rho * e     !! + (ω·θ·ρ0·e)

    end subroutine jwl_eos

    !!! Water EOS/Mie-Gruneisen
    !!! Type = 6
    elemental subroutine mie_gruneisen_eos_of_water(rho, e, p, c)
        real(8), intent(in)  :: rho, e
        real(8), intent(inout) :: p, c
        real(8) :: ratio, mu
        real(8), parameter :: rho0 = 1000, c0 = 1480, bulk = 2.2e9
        real(8), parameter :: gamma0 = 0.5, a = 0
        real(8), parameter :: s(3) = [2.56, 1.986, 1.2268]

        ratio = rho/rho0
        mu = ratio - 1

        if ( mu > 0 ) then
            p = rho0*c0**2*mu*(1+(1-gamma0/2)*mu-a/2*mu**2)       &
              / (1-(s(1)-1)*mu-s(2)*mu**2/ratio-s(3)*mu**3/ratio**2)**2 &
              + (gamma0+a*mu) * e
        else
            p = rho0*c0**2*mu + (gamma0+a*mu)*e
        end if

        c = sqrt(bulk / rho)

        ! if ( p < 0 ) p = 0

    end subroutine mie_gruneisen_eos_of_water

    !!! Water EOS/Polynomial
    !!! Type = 7
    elemental subroutine water_polynomial_eos(rho, e, p)
        real(8), intent(in)  :: rho, e
        real(8), intent(inout) :: p
        real(8) :: ratio, mu
        real(8), parameter :: rho0 = 1e3
        real(8), parameter :: a(3) = [2.19e9, 9.224e9, 8.767e9]
        real(8), parameter :: b(0:1) = [0.4934, 1.3937]

        ratio = rho/rho0
        mu = ratio - 1

        associate(poly => [1._8, mu, mu**2])
        if ( mu > 0 ) then
            p = dot_product(a, poly*mu) &
              + dot_product([b, b(1)], poly) * rho0 * e
        else
            p = a(1) * mu + dot_product(b, poly(1:2)) * rho0 * e
        end if
        end associate

    end subroutine water_polynomial_eos

    !!! Artificial equation of state for oil
    !!! Type = 8
    elemental subroutine oil_eos(rho, p)
        real(8), intent(in)  :: rho
        real(8), intent(inout) :: p
        real(8), parameter :: c = 50
        real(8), parameter :: rho0 = 917

        p = c**2 * (rho - rho0)

    end subroutine oil_eos
    
    !!! PETN gas EOS/Jones-Wilkins-Lee
    !!! Type = 9
    elemental subroutine jwl_eos_of_PETN(rho, e, p)
        real(8), intent(in)  :: rho, e
        real(8), intent(inout) :: p
        real(8) :: ratio
        real(8), parameter :: rho0 = 1770
        real(8), parameter :: A  = 6.17e11, B  = 16.9e9
        real(8), parameter :: R1 = 4.4,     R2 = 1.2
        real(8), parameter :: omega = 0.25

        ratio = rho/rho0

        p = A*(1-omega*ratio/R1)*exp(-R1/ratio) &
          + B*(1-omega*ratio/R2)*exp(-R2/ratio) &
          + omega * rho * e     !! + (ω·θ·ρ0·e)

    end subroutine jwl_eos_of_PETN

    !!! Solid EOS/Mie-Gruneisen For Armco Iron
    !!! Type = 101
    elemental subroutine mie_gruneisen_eos_of_armcoIron(rho, e, p)
        real(8), intent(in)  :: rho, e
        real(8), intent(inout) :: p
        real(8) :: mu
        real(8) :: p_H
        real(8), parameter :: rho0 = 7850
        real(8), parameter :: Gamma = 1.81

        mu = rho/rho0 - 1
        p_H = hugoniot_curve(mu)

        p = (1 - 0.5 * Gamma * mu) * p_H &
          + Gamma * rho * e

        !   if ( p < 0 ) p = 0

    contains
        elemental function hugoniot_curve(eta) result(re)
            real(8), intent(in) :: eta
            real(8) :: re
            real(8), parameter :: Cs = 3630
            real(8), parameter :: Ss = 1.8
            real(8) :: para(3)

            para(1) = rho0 * Cs**2
            para(2) = para(1) * (1 + 2 * (Ss - 1))
            para(3) = para(1) * (2 * (Ss - 1) + 3 * (Ss - 1)**2)

            if ( eta > 0 ) then
                re = sum(para*[eta, eta**2, eta**3])
            else
                re = para(1) * eta
            end if

        end function hugoniot_curve
    end subroutine mie_gruneisen_eos_of_armcoIron

    !!! Solid EoS
    !!! Type = 102
    elemental subroutine arti_eos_of_102(rho, p)
        real(8), intent(in)  :: rho
        real(8), intent(inout) :: p
        real(8), parameter :: rho0 = 1100, K = 2e7, cc = K / rho0

        p = cc * (rho - rho0)

    end subroutine arti_eos_of_102

    !!! Solid EoS
    !!! Type = 103
    elemental subroutine arti_eos_of_103(rho, p)
        real(8), intent(in)  :: rho
        real(8), intent(inout) :: p
        real(8), parameter :: rho0 = 2500, K = 0.33e6, cc = K / rho0

        p = cc * (rho - rho0)

    end subroutine arti_eos_of_103

    !!! Aluminium(Al5086) EoS
    !!! Type = 104
    elemental subroutine arti_eos_of_Aluminium(rho, p)
    real(8), intent(in)  :: rho
    real(8), intent(inout) :: p
    real(8), parameter :: rho0 = 2700, K = 7.235e10, cc = K / rho0

    p = cc * (rho - rho0)

end subroutine arti_eos_of_Aluminium

end module eos_m