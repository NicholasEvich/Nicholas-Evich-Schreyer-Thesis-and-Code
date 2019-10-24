function [dPdz_A] = dP_Acceleration(d, G, h_fg, qflux, v_fg, w)
% Function for calculating pressure DROP due to accelerational effects for 
% a straight pipe section of length dz
% 
% See equation __ on page __ of __
% 
%   d =  diameter (m)
%   G = mass flux (kg/s*m^2, mass flow rate divided by cross sectional area)
%   h_fg = latent heat of vaporization (J/kg)
%   qflux = heat flux (W/m^2)
%   v_fg = average specific volume of two-phase fluid (m^3/kg)
%   w = mass flow rate (kg/s)

dPdz_A = G * G * v_fg * (qflux * pi * d) / (w * h_fg);
end

% Consider using dx/dz (using x0 and x1 and z_inc) instead