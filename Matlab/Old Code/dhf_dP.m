function [diff] = dhf_dP(P0)
% Derivative of specific enthalpy (liquid water) with respect to pressure 
% (absolute). The function is the analytical derivative of a power
% function of saturated water data from NIST (____). It is not very
% accurate and should only be used as a placeholder for this function until
% I have better data or a better fit.

diff = 121543/P0;

end

