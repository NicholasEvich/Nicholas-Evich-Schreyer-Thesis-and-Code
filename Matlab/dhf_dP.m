function [diff] = dhf_dP(P0)
% Derivative of specific enthalpy (liquid water) with respect to pressure 
% (absolute). The function is the analytical derivative of a power
% function of saturated water data from NIST (____). It is not very
% accurate and should only be used as a placeholder for this function until
% I have better data or a better fit.
%   P0 is input in Pa, needs to be converted to MPa for this equation
P = P0/1E6;
diff = 1000*709.53*0.3525*P^(0.3525 - 1); % 1000 to go from kJ/kg to J/kg
end

