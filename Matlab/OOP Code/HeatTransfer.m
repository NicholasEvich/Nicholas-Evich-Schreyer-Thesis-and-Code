classdef HeatTransfer < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function obj = HeatTransfer()

        end
        
        % Break this function down into smaller ones
        function [h] = h(obj,d,w,qflux)
           if obj.x0 > 0.8 
              % flag error 
           end
           
           Re = obj.Re0;
           
           if Re < 2000
               % flag for laminar flow
           elseif Re > 5e6
               % flag for too turbulent flow
           end
           
           k_f = obj.K_F;
           pr = obj.PR;
           rho_f = obj.RHO_F;
           A = 0.25*pi*d*d;
           
           % should I assign these to the water class? is it even worth it?
           darcyf = (0.79*log(Re)-1.64)^-2;
           hsp = ((k_f/d)*(darcyf/8)*(Re - 1000)*pr)/(1 + 12.7*(darcyf/8)^(0.5)*(pr^(2/3) - 1));
           froude = (w/(A*rho_f))^(0.5)/(9.81*d);
           
           if froude >= 0.04
               ffr = 1;
           else
               ffr = 2.63*(froude)^(2/3);
           end
           
           % How inefficient is it to do this calculation twice?
           h1 = obj.h1(hsp,ffr,qflux,A,w);
           h2 = obj.h2(hsp,ffr,qflux,A,w);
               
           if h1 > h2
               h = h1;
           else
               h = h2;
           end       
        end
        
        function [h] = h1(obj, hsp, ffr, qflux, A ,w)
            rho_v = obj.RHO_V;
            rho_f = obj.RHO_F;
            x = obj.x0;
            h_fv = obj.H_FV;
            g_sf = obj.G_SF;
            
            h = hsp*(0.6683*(rho_f/rho_v)^0.1*x^0.16*(1 - x)^0.64*ffr + ...
                1058*(qflux*A/w/h_fv)^0.7*(1-x)^0.8*g_sf);
        end
        
        function [h] = h2(obj, hsp, ffr, qflux, A ,w)
            rho_v = obj.RHO_V;
            rho_f = obj.RHO_F;
            x = obj.x0;
            h_fv = obj.H_FV;
            g_sf = obj.G_SF;  
            
            h = hsp*(1.136*(rho_f/rho_v)^0.45*x^0.72*(1 - x)^0.08*ffr + ...
                667.2*(qflux*A/w/h_fv)^0.7*(1-x)^0.8*g_sf);
        end
    end
end

