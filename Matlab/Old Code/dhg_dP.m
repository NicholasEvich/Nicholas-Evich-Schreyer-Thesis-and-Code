function [diff] = dhg_dP(P0)
% Derivative of specific enthalpy (water vapor) with respect to pressure 
% (absolute). The function is the analytical derivative of a polyfit
% function of saturated water data from NIST (____). It is most accurate
% for pressures in excess of 50 kP (0.05 MPa)

diff = 0.0122*(2E6)*(P0)^(0.0122 - 1);

end

