classdef SaturatedWater < handle
    %UNTITLED2 Summary of this class goes here
    %   What if I only had one of these objects per pipe?
    
    % In reality, these aren't constant and may be changed later
    properties (Constant)
        % Saturation temperature of water (K)
        T_SAT = 373.15;
        % Latent heat of vaporization, valid at 373.15 K (J/kg)
        H_FV = 2256E03;
        % Density of saturated water at 373.15 K (kg/m^3)
        RHO_F = 958.05;
        V_F = 0.0010438;
        % Density of saturated steam at 373.15 K (kg/m^3)
        RHO_V = 0.590;
        V_G = 1.6949;
        % Specific heat of liquid water at constant pressure (J/kg*K)
        C_PF = 4181.3; %use NIST to make sure this is correct when I know the temperature
        % Thermal conductivity of water of saturdated water at 373.15 K (W/mK)
        K_F = 0.67703;
        % Dynamic viscosity of saturdated water at 373.15 K (Pa*s)
        MU_F = 0.0002814;
        % Dynamic viscosity of saturated steam at 373.15 K (Pa*s)
        MU_V = 11.97E-06;
        % Prandlt number for saturated water at 373.15 K (Dimensionless)
        PR = 1.76;
        % Difference between specific volumes of saturated water/steam at 373.15 K (m^3/kg)
        % V_DIFF = V_G - V_F;
        % Surface tension of saturated water at 373.15 K (N/m)
        SIGMA = 58.9e-03;
        % Multiplication factor for calculating heat transfer coefficient for water
        G_SF = 1;
    end
    
    properties
        T0 % current temperature
        P0 % current presssure v
        x0
        alpha0
%         T0 % current temperature
%         T_sat % This remains constant
%         H_fv % Latent ehat of vaporization
%         rho_f % density of fluid
%         rho_v % density of vapor
%         C_pf % specific heat of fluid phase
%         K_f % thermal conductivity of water at T0
%         mu_f % dynamic viscosity of water
%         mu_v % dynamic viscosity of vapor
%         PR % prandlt number
%         v_diff % difference in specific volumes
%         sigma % surface tension of water
%         G_sf % multiplication factor for calculating heat transfer coefficient
%         P0 % current pressure
    end
    
    methods
        function obj = SaturatedWater(T0, P0, x0)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.T0 = T0;
            obj.P0 = P0;
            obj.x0 = x0;
            
            rho_v = obj.RHO_V;
            rho_f = obj.RHO_F;
            
            obj.alpha0 = 1/(1 + rho_v/rho_f) * (1 - x0);

        end
        
        function [] = updatePressure(obj,drop)
            obj.P0 = obj.P0 - drop;
            if obj.P0 < 20000
                % Find a good way to flag what design this was and why its
                % a problem. In the future, maybe don't throw an error, but
                % discard this given object so that the rest of the
                % analysis can continue on without it
                error('Absolute pressure too low in pipe design.');
            end
        end
        
        % Need heat flux, diameter, dz, massflow rate for this one
        % will need to put this in another class or find another workaround
        % solution: the quality member belongs to this class, but the
        % function for updating it belongs to a derived class
        function [] = updateQuality(obj)
            obj.x0 = obj.x0;
        end
        
        % 12/3/19: This function is for calculating temperature changes for
        % a subcooled liquid. As I add more classes that take care of the
        % geometry, it may make more sense to NOT include qflux, d, dx, and
        % w in this function call. Also, C_PF may be unfavorable in
        % comparison to C_VF (constant volume)
        % In addition, I could probably include some branching
        % optimizations here. If used properly, I wouldn't even have to
        % include the branches. But, that would require trusting whoever is
        % developing that particular application
        % Should I be returning the temperature?
        function [] = updateTemp(obj, qflux, d, dz, w)
           if obj.T0 < obj.T_SAT
              deltaT = (qflux*d*pi*dz)/(obj.CP_F * w);
              obj.T0 = obj.T0 + deltaT;
           end
           
           if obj.T0 >= obj.T_SAT && obj.x0 < 1
              obj.T0 = obj.T_SAT; 
           end
           
%            % If I am going to do it this way, I need to make sure that
%            % temperature isnt used for anything before it updates again
%            if obj.T0 >= obj.T_SAT && obj.x0 < 1
%               obj.T0 = obj.T_SAT;
%            elseif obj.T0 < obj.T_SAT
%               deltaT = (qflux*d*pi*dz)/(obj.CP_F * w);
%               obj.T0 = obj.T0 + deltaT;
%            end
        end

% The four functions below are for calculating the approximate values of
% the derivatives of specific volume and enthalpy, respectively, with
% respect to absolute pressure. These functions were interpolated using
% various methods on data from "Thermophysical Properties of Fluid Systems"
% at NIST. Specifically, the "Saturation properties — temperature
% increments" option was selected. These functions are not perfect, and
% significant errors can exist at very high or very low pressures. 
% ----------- See https://webbook.nist.gov/chemistry/fluid/ --------------
% These functions require the usage of MPa, so either utilize them such
% that only pressures in MPa are used, or make the appropriate conversions
% within the functions themselves.

% For future versions, 
%       get more consistent with units, get better functions (more
%       accurate), define constraints on accuracy of functions
        
        function outputArg = dvf_dP(obj)
            P = obj.P0*1e-06;
            outputArg = -18E-7*(P)^5 + 20E-6*(P)^4 -12E-5*(P)^3 + 3E-4*(P)^2 - 0.0004*P + 0.0002;
        end
        
        function outputArg = dvg_dP(obj)
            P = obj.P0*1e-06;
            outputArg = -0.96*(0.1909)*P^(-1.96);
        end
        
        function outputArg = dhf_dP(obj)
            P = obj.P0*1e-06;
            outputArg = 1000*709.53*0.3525*P^(0.3525 - 1); % 1000 to go from kJ/kg to J/kg
        end
        
        function outputArg = dhg_dP(obj)
            P = obj.P0*1e-06;
            outputArg = 1000*30.7/P; % 1000 to go from kJ/kg to J/kg
        end
        
        % Again, think about ways that I can avoid the geometry of the
        % situation here. It is entirely possible that this is the best way
        % to do it, but I don't know
        function [] = quality_alpha(obj, d, qflux, dz, w)
            C_pf = obj.C_PF;
            T_sat = obj.T_SAT;
            T = obj.T0; % current temperature
            h_fv = obj.H_FV;
            
            obj.x0 = (-C_pf*(T_sat - T))/h_fv + (pi*d*qflux*dz)/(w*h_fv) + obj.x0;
            
            % Pick one of these, both of them are completely unnecessary
%             if obj.x0 > 1
%                 obj.x0 = 1;
%             end
            
            if obj.x0 > 0.8
                error('Flow quality too high to compute accurate heat transfer coefficients');
            end
            
            rho_v = obj.RHO_V;
            rho_f = obj.RHO_F;
            
            obj.alpha0 = 1/(1 + rho_v/rho_f) * (1 - obj.x0);
            
        end
    end
end

