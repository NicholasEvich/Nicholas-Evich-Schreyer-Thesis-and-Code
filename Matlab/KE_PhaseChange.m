function [factor] = KE_PhaseChange(G, h_fg, v_f, v_g, v_fg, xe)
% The output of this function represents the ratio of kinetic energy
% changes resulting from phase change due to the latent heat exchange from
% evaporation or condensation
%  See equation __ on page __ of __
%     G = mass flux (kg/s*m^2, divide mass flow rate by cross sectional area)
%     h_fg = latent heat of vaporization (J/kg)
%     v_f = specific volume of liquid (m^3/kg)
%     v_g = specific volume of vaport (m^3/kg)
%     v_fg = average specific volume of fluid mixture (m^3/kg)
%     xe = thermodynamic equilibrium quality (is assumed to be equal to the
%     normal flow quality when the two fluid phases are in thermal
%     equilibrium

factor = 1 + G * G * v_fg / h_fg * (xe * v_g + (1 - xe) * v_f);
end

