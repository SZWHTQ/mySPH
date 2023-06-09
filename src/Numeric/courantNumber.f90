module cour_num_m
    use ctrl_dict, only: Project, Config
    implicit none

contains

    elemental function courant_num(hsml, div_v, c) result(xi)
        real(8), intent(in)  :: hsml, div_v, c
        real(8) :: xi
        real(8) :: alpha, beta


        select case(Project%nick)
            ! case ("dam_break")
            !     alpha = 1
            !     beta  = 1
            !     psi   = 0.1
            ! case("shear_cavity","shock_tube")
            !     alpha = 1
            !     beta  = 1
            !     psi   = 0.1
            case("tnt_bar", "tnt_cylinder", "undex_cylinder", "undex_chamber")
                alpha = 1
                beta  = 10
            case("taylor_rod")
                alpha = 0.5
                beta  = 0.5
            case("can_beam", "can_beam_3d", "can_beam_f")
                alpha = 0.5
                beta  = 0.5
            case("beam_oil")
                alpha = 0.5
                beta  = 0.5
            case("UNDEX")
                alpha = 1
                beta  = 1
            case("water_impact")
                alpha = 1
                beta  = 1
                ! if ( P(i)%Type > 100 ) then
                !     alpha = 1
                !     beta  = 1
                ! end if
            case("undex_plate")
                alpha = 1
                beta  = 1
                ! if ( P(i)%Type > 100 ) then
                !     alpha = 2
                !     beta  = 2
                ! end if
            case default
                alpha = 1
                beta  = 1
            end select

        xi = Config%delta_t &
            * ( hsml * div_v + c + 1.2*(alpha*c+beta*hsml*abs(div_v)) ) &
            / hsml

    end function courant_num

end module cour_num_m