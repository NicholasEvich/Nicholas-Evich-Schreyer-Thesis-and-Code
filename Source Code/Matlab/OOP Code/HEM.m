classdef HEM < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function obj = HEM()

        end

        function total_drop = dPdz(obj, theta, d, qflux, w)
            G = w/(0.25*pi*d*d);
            
            if G > 1e3
               % flag error but do not exit code
               % how do I flag this error within the HEM class? The error
               % needs to be flagged in the actual Pipe object
            end
               
            KE = obj.KE_PhaseChange(G);
            Grav = obj.dP_gravity(theta);
            Fr = obj.dP_friction(d, G);
            A = obj.dP_acceleration(d, G, qflux, w);
            C = obj.compressibility(G);
            F = obj.flashing(G);
            total_drop = (KE * Fr + A + Grav) / (KE);
        end
        
        % This works, could probably use some optimizations though
        function drop = dP_gravity(obj, theta)
            v_fg = obj.V_G - obj.V_F;
            v_f = obj.V_F;
            x = obj.x0;
            drop = 9.81 * sind(theta) / (v_f + x * v_fg);
        end
        
        function drop = dP_friction(obj, d, G)
            v_fg = obj.V_G - obj.V_F;
            v_f = obj.V_F;
            x = obj.x0;
            mu_f = obj.MU_F;
            mu_v = obj.MU_V;
            
            % replace name placeholder once I know what this physically
            % represents
            placeholder = G*d/mu_f;
            
            if placeholder < 2300
                c = 16;
                N = 1;
            elseif placeholder > 20000
                c = 0.046;
                N = 0.20;
            else
                c = 0.079;
                N = 0.25;
            end
            
            ffo = c/(placeholder)^N;
            mu_bar = x*mu_v + (1-x)*mu_f;
            
            fricFac = ffo*(mu_bar/mu_f)^N;
            
            drop = (2 / d) * fricFac * v_f * G * G * (1 + x * v_fg / v_f);
        end
        
        function drop = dP_acceleration(obj, d, G, qflux, w)
            v_fg = obj.V_G - obj.V_F;
            h_fg = obj.H_FV;
            drop = G * G * v_fg * (qflux * pi * d) / (w * h_fg);
        end
        
        function K = KE_PhaseChange(obj, G)
            v_f = obj.V_F;
            v_g = obj.V_G;
            v_fg = v_g - v_f;
            h_fg = obj.H_FV;
            % Be explicit on the difference and usages of x and xe
            xe = obj.x0;
            K = 1 + G * G * v_fg / h_fg * (xe * v_g + (1 - xe) * v_f);
        end
        
        function K = compressibility(obj, G)
            xe = obj.x0;
            K = G*G*(xe * obj.dvg_dP() + (1 - xe) * obj.dvf_dP());
            % G*G*(xe*dvg_dP(p_abs) + (1 - xe)*dvf_dP(p_abs));
        end
        
        function K = flashing(obj, G)
            v_fg = obj.V_G - obj.V_F;
            h_fg = obj.H_FV;
            xe = obj.x0;
            K = (G * G* v_fg / h_fg) * (xe * obj.dhg_dP() + (1 - xe) * obj.dhf_dP());
            % (G*G*v_fg/h_fg)*(xe*dhg_dP(p_abs) + (1 - xe)*dhf_dP(p_abs));
        end
    end
end

