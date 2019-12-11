classdef HEM < SaturatedWater
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    % Does this class even require any properties
    properties
        
    end
    
    methods
        function obj = HEM(T0, P0, x0)
            % T0, P0, and x0 are the entrance temperature, pressure, and
            % flow quality for the fluid channel or object being analyzed.
            % The idea here is that this class (and the SaturatedWater
            % class that it inherits from) can function completely
            % independently of the geometry of the flow scenario. Any
            % geometric factors that need to be considered can be handled
            % at the application level or in a separate channel/pipe class.
            % Therefore, some of the formulas and functions used here may
            % omit any terms that depend on the geometry.
            obj@SaturatedWater(T0, P0, x0);
        end
        
        % This class and its functions only calculate pressure DROPS. The
        % value of the absolute pressure is handled by other classes
        % This function should know nothing of dz. That is not the job of
        % this class
        
        % This function call can definitely be optimized: have a function
        % for calculating the friction factor, etc.
        function total_drop = dPdz(obj, theta, G, d, fricFac, qflux, w)
            KE = obj.KE_PhaseChange(G);
            G = obj.dP_gravity(theta);
            Fr = obj.dP_friction(d, fricFac, G);
            A = obj.dP_acceleration(d, G, qflux, w);
            C = obj.compressibility(G);
            F = obj.flashing(G);
            total_drop = (KE * Fr + A + G) / (KE + C - F);
        end
        
        % This works, could probably use some optimizations though
        function drop = dP_gravity(obj, theta)
            v_fg = obj.V_G - obj.V_F;
            v_f = obj.V_F;
            x = obj.x0;
            drop = 9.81 * sind(theta) / (v_f + x * v_fg);
        end
        
        function drop = dP_friction(obj, d, fricFac, G)
            v_fg = obj.V_G - obj.V_F;
            v_f = obj.V_F;
            x = obj.x0;
            drop = (2 / d) * fricFac * v_f * G * G * (1 + x * v_fg / v_f);
        end
        
        % Is it really necessary to pass d, G, and w together? I could
        % probably get away with just two of those ut I would need to have
        % something like pi/4 already calculated as a constant data member
        % or something like that
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
        
        % I am not exactly sure if I am using the obj.superclass_method()
        % syntax correctly...
        function K = compressibility(obj, G)
            xe = obj.x0;
            K = G*G*(xe * obj.dvg_dP() + (1 - xe) * obj.dvf_dP());
        end
        
        function K = flashing(obj, G)
            v_fg = obj.V_G - obj.V_F;
            h_fg = obj.H_FV;
            xe = obj.x0;
            K = (G * G* v_fg / h_fg) * (xe * obj.dhg_dP() + (1 - xe) * obj.dhf_dP());
        end
        
        % Next functions: (really think about whether these functions
        % belong in the HEM class or the saturated water class)
        %   friction factor for pressure drop equations 
    end
end

