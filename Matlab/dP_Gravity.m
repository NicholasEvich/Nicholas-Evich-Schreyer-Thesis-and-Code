function [dPdz_G] = dP_Gravity(theta, v_f, v_fg, x)
% Function for calculating pressure DROP due to gravitational effects for a
% straight pipe section of length dz
% 
% See equation __ on page __ of __
% 
%   theta = incline of pipe (degrees)
%   v_f = specific volume of liquid phase (m^3/kg)
%   v_fg = average specific volume of two-phase fluid (m^3/kg)
%   x = vapor quality
g = 9.81;

dPdz_G = g * sind(theta) / (v_f + x * v_fg);
end

