function [diff] = dvg_dP(P0)
% Derivative of specific volume (water vapor) with respect to pressure 
% (absolute). The function is the analytical derivative of a polyfit
% function of saturated water data from NIST (____). It is most accurate
% for pressures in excess of 50 kP (0.05 MPa)
%   P0 is input in Pa, needs to be converted to MPa for this equation
P = P0/1E6;
diff = -0.96*(0.1909)*P^(-1.96);
end

