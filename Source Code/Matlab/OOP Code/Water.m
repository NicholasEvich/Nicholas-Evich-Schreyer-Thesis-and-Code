 classdef Water < HEM & HeatTransfer

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
        u0 % single-phase velocity
        Re0 % Reynolds number
        darcyf0 % darcy friction factor
        hsp0 % single phase convection coefficient
        froude0 % Froude number
        phase % string or enumeration for the various states
    end
    
    methods
        function obj = Water(T0, P0, x0)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj@HEM();
            obj@HeatTransfer();
            
            obj.T0 = T0;
            obj.P0 = P0;
            obj.x0 = x0; % this should be calculated separately
            
            rho_v = obj.RHO_V;
            rho_f = obj.RHO_F;
            
            obj.alpha0 = 1/(1 + rho_v/rho_f * (1 - x0)/x0);

        end
        
        function [] = updateRe(obj, d, w)
            rho_f = obj.RHO_F;
            mu_f = obj.MU_F;
            
            A = 0.25*pi*d*d;
            u = w/(rho_f*A);
            obj.u0 = u;
            obj.Re0 = rho_f*u*d/mu_f;
        end
        
        function [] = updatePressure(obj, dz, drop)
            obj.P0 = obj.P0 - dz*drop;
        end

        function [] = updateTemp(obj, qflux, d, dz, w)
           if obj.T0 < obj.T_SAT
              deltaT = (qflux*d*pi*dz)/(obj.CP_F * w);
              obj.T0 = obj.T0 + deltaT;
           end
           
           if obj.T0 >= obj.T_SAT && obj.x0 < 1
              obj.T0 = obj.T_SAT; 
           end
           
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
        
% Add ways to ensure these functions only at upon saturated water,
% not subcooled versions
        function outputArg = dvf_dP(obj)
            P = obj.P0*1e-06;
            outputArg = -18E-7*(P)^5 + 20E-6*(P)^4 -12E-5*(P)^3 + 3E-4*(P)^2 - 0.0004*P + 0.0002;
        end
        
        function outputArg = dvg_dP(obj)
            P = obj.P0*1e-06;
            outputArg = -0.96*(0.1909)*P^(-1.96);
        end
        
        function outputArg = dhf_dP(obj)
            outputArg = 121543/obj.P0;
        end
        
        function outputArg = dhg_dP(obj)
            outputArg = 0.0122*(2E6)*(obj.P0)^(0.0122 - 1);
        end
        
        function [] = updateQuality(obj, d, qflux, dz, w)
            C_pf = obj.C_PF;
            T_sat = obj.T_SAT;
            T = obj.T0; % current temperature
            h_fv = obj.H_FV;
            
            obj.x0 = (-C_pf*(T_sat - T))/h_fv + (pi*d*qflux*dz)/(w*h_fv) + obj.x0;
            
            rho_v = obj.RHO_V;
            rho_f = obj.RHO_F;
            
            obj.alpha0 = 1/(1 + rho_v/rho_f) * (1 - obj.x0);
            
            % There is probably a better way to do this while
            % simultaneously updating the fluid phase
            if obj.x0 >= 1
               obj.x0 = 1; obj.alpha0 = 1; 
            end
            
        end
    end
end

