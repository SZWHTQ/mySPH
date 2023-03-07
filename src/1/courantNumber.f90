module cour_num_m
    use ctrl_dict, only: delta_t
    use parse_toml_m, only: nick
    implicit none

contains

    elemental function courant_num(hsml, div_v, c) result(xi)
        real(8), intent(in)  :: hsml, div_v, c
        real(8) :: xi
        real(8) :: alpha, beta

        select case(nick)
        case("dam_break")
            alpha = 1
            beta  = 1
        case("shear_cavity","shock_tube")
            alpha = 1
            beta  = 1
        case("tnt_bar", "tnt_cylinder", "undex_cylinder", "undex_chamber")
            alpha = 1
            beta  = 10
        case default
            alpha = 1
            beta  = 1
        end select

        xi = delta_t &
        * ( hsml * div_v + c + 1.2*(alpha*c+beta*hsml*abs(div_v)) ) &
        / hsml

    end function courant_num

end module cour_num_m