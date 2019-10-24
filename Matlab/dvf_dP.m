function [diff] = dvf_dP(P0)
% Derivative of specific volume (liquid water) with respect to pressure 
% (absolute). The function is the analytical derivative of a polyfit
% function of saturated water data from NIST (____). It is most accurate
% for pressures in excess of 50 kP (0.05 MPa)
%   P0 is input in Pa, needs to be converted to MPa for this equation
P = P0/1E6;
diff = -18E-7*(P)^5 + 20E-6*(P)^4 -12E-5*(P)^3 + 3E-4*(P)^2 - 0.0004*P + 0.0002;
end

