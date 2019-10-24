function [diff] = dhg_dP(P0)
% Derivative of specific enthalpy (water vapor) with respect to pressure 
% (absolute). The function is the analytical derivative of a polyfit
% function of saturated water data from NIST (____). It is most accurate
% for pressures in excess of 50 kP (0.05 MPa)
%   P0 is input in Pa, needs to be converted to MPa for this equation
P = P0/1E6;
diff = 1000*30.7/P; % 1000 to go from kJ/kg to J/kg
end

