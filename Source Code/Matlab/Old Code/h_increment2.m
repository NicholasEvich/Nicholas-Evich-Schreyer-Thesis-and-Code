function [h_inc2] = h_increment2(hsp, rho_f, rho_v, x, ffr, qflux, A, H_fv ,w, G_sf)
% This function is used for calculating the heat transfer coefficient for a
% straight pipe. In practice, it will be used on a very small segment of a
% pipe many times in order to calculate an overall heat transfer
% coefficient. The value that will actually be used will be the higher of
% the return values for h_increment1 and h_increment 2.
% -----------------------------------------------------------------------
% See equation _ on page _ of _
% -----------------------------------------------------------------------
%     hsp = Single phase convection coefficient
%     rho_f = Density of fluid phase
%     rho_v = Density of vapor phase
%     x = instantaneous flow quality
%     ffr = stratification parameter, a function of the Froude Number
%     qflux = Average heat flux for a given pipe design
%     A = Instantaneous cross sectional area of the pipe
%     H_fv = Latent heat of vaporization
%     w = Mass flow rate
%     G_sf = 

h_inc2 = hsp*(1.136*(rho_f/rho_v)^0.45*x^0.72*(1 - x)^0.08*ffr + ...
    667.2*(qflux*A/w/H_fv)^0.7*(1-x)^0.8*G_sf);
end

