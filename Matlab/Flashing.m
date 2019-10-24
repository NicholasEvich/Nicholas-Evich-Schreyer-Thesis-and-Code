function [K_flashing] = Flashing(G, h_fg, p_abs, v_fg, xe)
% The output of this function represents the contribution to overall
% pressure drop due to compressibility effects. This effect is amplified as
% more fluid is converted to vapor.
% See equation __ on page __ of __
%    G = mass flux (kg/s*m^2, divide mass flow rate by cross sectional area)
%    p_abs = absolute pressure
%    xe = thermodynamic equilibrium quality (is assumed to be equal to the
%    normal flow quality when the two fluid phases are in thermal
%    equilibrium
K_flashing = (G*G*v_fg/h_fg)*(xe*dhg_dP(p_abs) + (1 - xe)*dhf_dP(p_abs));
% Replace h_fg and v_fg with calls to those variables as a function of
% pressure
end

