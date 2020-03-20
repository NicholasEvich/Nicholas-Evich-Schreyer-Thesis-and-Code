classdef Pipe < handle
    properties (Constant)
        d_min = 0.005
        d_max = 0.032
    end
    
    properties
        N % number of discrete channel segments to evaluate
        L % total length of channel (m)
        d0 % inlet diameter, may not be used
        % phi
        % theta
        w % (constant) mass flow rate of fluid through channel
        dz % incremental segment length in meters, given by L/N
        diams % length-N array of incremental diameters based on diameter function
        thetas % length-N array of incremental inclines based on centerline function
        SA % total surface area of channel, calculated through integration of perimeter function (derived from diameter function)
        fluidProps % object of class <fluid_name>, used to evaluate pressure drops and heat transfer coefficients
        
        Q % total heat input (this will be replaced by a heat flux function later on
        qfluxes % array of heat fluxes
    end
    
    methods
        % How can I make the declarations below simpler (like JavaScript)?
        % And, how can I declare constant properties in the constructor?
        % I'd like to make the angles constant, but I still need the
        % freedom to actually initialize them here
        function obj = Pipe(L, N, w, r, centerline, qflux) 
            diams = zeros(1,N);
            thetas = zeros(1,N);
            fluxes = zeros(1,N);
            
            % evaluating all incremental diameters and inclines
            % this loop also calculates the surface area of the channel
            % through numeric integration
            
            dz = L/N;
            
            for i = 1:N
                z = i*dz; % which is more efficient, this or z += dz;?
                diams(i) = 2*r(z);
                % for now, make all theta functions = 0
                % SA = SA + dz*2*pi*d(z); % verify my formula for this
                slope = (centerline(z + dz) - centerline(z - dz))/(2*dz); % three-point centered difference formula
                thetas(i) = atand(slope);
                
                fluxes(i) = qflux(z);
            end  
            
            P = @(z) 2*pi*r(z);
            
            SA = comp_trapz(P,0,L,N);
            
            % storing all incremental diameters, inclines, and total
            % surface area
            obj.L = L;
            obj.N = N;
            obj.dz = L/N;
            obj.w = w;
            obj.diams = diams;
            obj.thetas = thetas;
            obj.qfluxes = fluxes;
            
            % obj.Q = Q;
            obj.SA = SA;
        end
        
        function [] = initFluidProps(obj,T_inlet, P_inlet, x_inlet)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            % perform more rigorous verification of input correctness in
            % the actual fluid (water class)- certain combinations of T, P,
            % and x clearly are not possible
            obj.fluidProps = Water(T_inlet, P_inlet, x_inlet);
        end
        
        % add h_overall to this
        function [deltaP, h_overall, x_out] = evalDesign(obj)
            fluid = obj.fluidProps;
            % the naming scheme and comments below are okay, maybe make it
            % less confusing though
            diams = obj.diams; %#ok<*PROP>
            thetas = obj.thetas; %#ok<*PROP>
            deltaPdz = 0;
            h_total = 0;
            w = obj.w; %#ok<*PROP>
            dz = obj.dz; %#ok<*PROP>
            % the qflux below will be replaced by a flux vector for each d
            qfluxes = obj.qfluxes; %#ok<*PROP>

            for i = 1:obj.N
               % pressure drops, quality changes, and heat transfer coef
               % do I store all this data in this object or in the water
               % object?
               d = diams(i);
               theta = thetas(i);
               qflux = qfluxes(i);

               fluid.updateRe(d,w);
               
               inc_drop = fluid.dPdz(theta, d, qflux, w);
               fluid.updatePressure(dz, inc_drop); % updates fluid object
               fluid.updateQuality(d, qflux, dz, w);
               deltaPdz = deltaPdz + inc_drop; % sums up total pressure drop
               
               h_inc = fluid.h(d,w,qflux);
               h_total = h_total + h_inc;
            end
            
            deltaP = dz*deltaPdz;
            
            h_overall = h_total/obj.N;
            
            % fluid.P0 = fluid.P0 - deltaP;
            if fluid.P0 < 0
               % flag an error/infeasibility 
            end
                      
            x_out = fluid.x0;
        end
    end
end

