function [K_compressibility] = Compressibility(G, p_abs, xe)
% The output of this function represents the contribution to overall
% pressure drop due to compressibility effects. This effect is amplified as
% more fluid is converted to vapor.
% See equation __ on page __ of __
%    G = mass flux (kg/s*m^2, divide mass flow rate by cross sectional area)
%    p_abs = absolute pressure
%    xe = thermodynamic equilibrium quality (is assumed to be equal to the
%    normal flow quality when the two fluid phases are in thermal
%    equilibrium
K_compressibility = G*G*(xe*dvg_dP(p_abs) + (1 - xe)*dvf_dP(p_abs));
end

