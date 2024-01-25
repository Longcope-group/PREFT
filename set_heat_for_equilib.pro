; given a tube structure with density and temperature set, set tube.heat
;  to balance radiative losses and thermal conduction, thus
;  thus puting the tube in equilibrium

pro set_heat_for_equilib, tube

dfc = calc_dfc_direct( tube )
calc_rad_loss, tube

tube.heat = ( tube.rad_loss - tube.p*dfc*tube.dl/tube.t/(tube.gam-1.0) ) > 0.0

return
end