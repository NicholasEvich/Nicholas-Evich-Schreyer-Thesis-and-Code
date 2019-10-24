function [dPdz_F] = dP_Friction(d, fricFac, G, v_f, v_fg, x)
% Function for calculating pressure DROP due to frictional effects for a
% straight pipe section of length dz
% 
% See equation __ on page __ of __
% 
%   d = hydraulic diameter (m)
%   fricFac = friction factor (related to wall shear stress)
%   G = mass flux (kg/s*m^2, mass flow rate divided by cross sectional area)
%   v_f = specific volume of liquid phase (m^3/kg)
%   v_fg = average specific volume of two-phase fluid (m^3/kg)
%   x = vapor quality

dPdz_F = (2 / d) * fricFac * v_f * G * G * (1 + x * v_fg / v_f);
end

