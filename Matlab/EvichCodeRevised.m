%

clc, clear;

% ---------------- Defining all constants-------------------------- %
% --------- ALL CONSTANTS BEGIN WITH A CAPITAL LETTER ------------- %
% Length of pipe (m)
L = 1;
% Gravitational constant (m/s^2)
G = 9.81;
% Inlet temperature of fluid (K)
T_inlet = 373.15;
% Saturation temperature of water (K)
T_sat = 373.15;
% Latent heat of vaporization, valid at 373.15 K (J/kg)
H_fv = 2256E03;
% Density of saturated water at 373.15 K (kg/m^3)
RHO_f = 958.05;
% Density of saturated steam at 373.15 K (kg/m^3)
RHO_v = 0.590;
% Specific heat of liquid water at constant pressure (J/kg*K)
C_pf = 4181.3; %use NIST to make sure this is correct when I know the temperature
% Thermal conductivity of water of saturdated water at 373.15 K (W/mK)
K_f = 0.67703;
% Dynamic viscosity of saturdated water at 373.15 K (Pa*s)
MU_f = 0.0002814;
% Dynamic viscosity of saturated steam at 373.15 K (Pa*s)
MU_v = 11.97E-06;
% Prandlt number for saturated water at 373.15 K (Dimensionless)
PR = 1.76;
% Average specific volume of saturated water/steam at 373.15 K (m^3/kg)
V_avg = 1 / RHO_v - 1 / RHO_f;
% Surface tension of saturated water at 373.15 K (N/m)
SIGMA = 58.9e-03;
% Multiplication factor for calculating heat transfer coefficient for water
G_sf = 1;
% Vertical incline of pipe
THETA = 0;

% --------------- Independent variables -------------------------
% Mass flow rate (kg/s)
% massflow = linspace(0.5,1.5,10);
massflow = 0.5;
% Inlet diameter satisfying confinement number (m)
% inlet_diameter = linspace(0.006,0.1, 10);
inlet_diameter = 0.01;
% Pipe angle in degrees
pipe_angle = linspace(0,45,100); 
% Heat inputs from 500 W to 5 MW in multiples of 10
% heat_input = 5*logspace(2,6,5);
heat_input = 500;
%
z = linspace(0, L, 250); % Consider making this an odd number so I get the middle section
% or use some other mathematical method to reduce error
z_increment = z(2) - z(1);

% ------------ Preallocation of memory for arrays ---------------

quality = zeros(length(massflow), length(pipe_angle), length(inlet_diameter), length(heat_input), length(z));
quality_h = zeros(length(massflow), length(pipe_angle), length(inlet_diameter), length(heat_input), length(z));
quality_p = zeros(length(massflow), length(pipe_angle), length(inlet_diameter), length(heat_input), length(z) - 1);
alpha = zeros(length(massflow), length(pipe_angle), length(inlet_diameter), length(heat_input), length(z));
% Order this one in a way that I can easily few the data when selecting the
% variable
diam_z = zeros(length(inlet_diameter),length(pipe_angle), length(z));
surface_area = zeros(length(inlet_diameter), length(pipe_angle));
cross_area = zeros(length(inlet_diameter), length(pipe_angle), length(z));
% h_inc = zeros(length(z));
h_inc = linspace(0, 0, 250); % Why is this necessary
h = zeros(length(massflow), length(pipe_angle), length(inlet_diameter), length(heat_input));
pressure_inc = linspace(0, 0, 249); % Why is this necessary
pressure_drop = zeros(length(massflow), length(pipe_angle), length(inlet_diameter), length(heat_input));

%
% ----------------- Calculating diameter as a function of z for every pipe
% angle ahead of time (this way it is not done every time there is a new
% heat input, inlet diameter, or massflow)
%       This loop also avoids any and all conditionals (saves time)

%Do = 0.1;
for i = 1:length(inlet_diameter)
    for j = 1:length(pipe_angle)
        
        surface_area(i,j) = pi*inlet_diameter(i)*L + pi*L^2*tand(pipe_angle(j))/2;
        
        for k = 1:length(z)/2
            diam_z(i,j,k) = inlet_diameter(i) + 2*z(k)*tand(pipe_angle(j));
            cross_area(i,j,k) = 0.25*pi*diam_z(i,j,k)*diam_z(i,j,k);
        end

        for k = length(z)/2:length(z)
            diam_z(i,j,k) = inlet_diameter(i) + 2*L*tand(pipe_angle(j)) - 2*z(k)*tand(pipe_angle(j));
            cross_area(i,j,k) = 0.25*pi*diam_z(i,j,k)*diam_z(i,j,k);
        end
    end
end

for i = 1:length(massflow)
    w = massflow(i);
    for j = 1:length(pipe_angle)
        phi = pipe_angle(j);
        for k = 1:length(inlet_diameter)
            d0 = inlet_diameter(k);
            for m = 1:length(heat_input)
                q_flux = heat_input(m)/surface_area(k,j); % I may be able 
                % change the order of these loops in order to avoid
                % initializing this variable too many times
                % Quality(i,j,k,m,1) is supposed to equal zero, so I dont
                % need to change anything since I already initialized it to
                % that
                for n = 2:length(z)
                    z1 = z(n); % For future versions, just use z(n) and z(n-1)
                    z0 = z(n-1);
                    D = diam_z(k,j,n);
                    A = cross_area(k,j,n); % Same as above
                    
                    %Xe(count1, count2, count3, count4, count5) = ...
                    %    (-Cp_f*(Tsat-Ti))/hfg + (pi*D*q*z(count5))/(W*hfg);
                    % In Nick's code, his quality calculations are the same
                    % (I think), yet he gets different values for each one
                    % But, my regular (utopia) quality calcs are "correct"
                    % Replace (z1 -z0) with a global z_increment constant
                    quality(i,j,k,m,n) = (-C_pf*(T_sat - T_inlet))/H_fv + ...
                        (pi*D*q_flux*(z1 - z0))/(w*H_fv) + quality(i,j,k,m,n-1);
                    quality_h(i,j,k,m,n) = quality(i,j,k,m,n); % Optimize this
                    % These are the "correct" values for both arrays
                    
                    if quality(i,j,k,m,n) > 1
                        quality(i,j,k,m,n) = 1;
                    end
                    
                    alpha(i,j,k,m,n) = 1/(1 + RHO_v/RHO_f)*...
                        (1 - quality(i,j,k,m,n));
                    
                    if quality_h(i,j,k,m,n) > 0.8
                       quality_h(i,j,k,m,n) = NaN; 
                    end
                    
                    % ---------------------------------------------------
                    % Everything after this can use further optimization
                    
                    u = w/(RHO_f.*A); % Calculating velocity of single phase (m/s)
                    Re = RHO_f*u*D/MU_f; % Reynolds number at specific conditions
               
                    if Re < 2000
                        error(['Flow cannot be laminar. Change input'...
                            'conditions, i.e. increase massflow or '...
                            'decrease diameter.'])
                    elseif Re > 5e6
                        error(['Flow is too turbulent. Change input'...
                            'conditions, i.e decrease massflow or '...
                            'increase diameter.'])
                    end
                    
                    darcyf = (0.79*log(Re)-1.64)^2; % Calculating darcy friction factor
                    % for smooth pipes
                    hsp = (K_f/D)*(darcyf/8*(Re - 1000)*PR)/...
                        (1 + 12.7*(darcyf/8)^0.5*(PR^(2/3)-1)); % Calculating
                    % single phase convection coefficient
                    
                    if (w/(A*RHO_f))^2/(G*D) >= 0.04 % Conditional statement to 
                        % determine FFr (what is FFr)
                        FFr = 1;
                    else % Figure out where this equation comes from and make s
                        % sure that the order of operations is working
                        % correctly
                        FFr = 2.63*((w/(A*RHO_f))^2/(G*D))^0.3;
                    end
                    
                    % Heat transfer coefficient increment eq 1 (check
                    % OoO's)
                    h_inc1 = hsp*(0.6683*(RHO_f/RHO_v)^0.1*...
                        quality_h(i,j,k,m,n)^0.16*...
                        (1-quality_h(i,j,k,m,n))^0.64...
                        *FFr+1058*(q_flux*A/w/H_fv)^0.7*(1-... % Whats up with all these divisions?
                        quality_h(i,j,k,m,n))^0.8*G_sf);
                    % Heat transfer coefficient increment eq 2
                    h_inc2 = hsp*(1.136*(RHO_f/RHO_v)^0.45*...
                        quality_h(i,j,k,m,n)^0.72*...
                        (1-quality_h(i,j,k,m,n))^0.08...
                        *FFr+667.2*(q_flux*A/w/H_fv)^0.7*(1-... % Whats up with all these divisions?
                        quality_h(i,j,k,m,n))^0.8*G_sf);
                    if h_inc1 >= h_inc2 
                        h_inc(n) = h_inc1;
                    else
                        h_inc(n) = h_inc2;
                    end
                    
                    % Why is this here?
                    if sqrt(SIGMA/(G*(RHO_f-RHO_v)))/D >= 0.5 % Setting
                        % confinement number limit
                        h(i,j,k,m) = NaN;
                    else
                        h(i,j,k,m) = mean(h_inc);
                    end
                    
                    %-----------------------------------------------------
                end
                
                for n = 1:(length(z) - 1)
                    % z1 = z(n); % For future versions, just use z(n) and z(n-1)
                    % z0 = z(n - 1);
                    D = diam_z(k,j,n); % This is the second time I am calculating these,
                    A = cross_area(k,j,n); % so precalculate them and store in an array
                    quality_p(i,j,k,m,n + 1) = (-C_pf*(T_sat - T_inlet))/H_fv + ...
                        (pi*D*q_flux*(z_increment))/(w*H_fv) + quality_p(i,j,k,m,n);
                    
                    % Ideally, remove these conditionals and use some kind
                    % of array to store c and n values
                    G = w/A; %Calculating mass flow rate per area(kg/m^2*s)
                    if G*D/ MU_f < 2300 %Determining c and n constants for two-phase friction factor
                        c = 16;
                        N = 1; % Capital N to avoid mixing up loop index (don't use n in the future)
                    elseif G*D/MU_f > 20000
                        c = 0.046;
                        N = 0.20;
                    else 
                        c = 0.079;
                        N = 0.25;
                    end
                    
                    ffo = c/(G*D/MU_f)^N; 
                    MU_bar = quality_p(i,j,k,m,n)*MU_v...
                        + (1-quality_p(i,j,k,m,n))*MU_f;
                    fTP = ffo*(MU_bar/MU_f)^N; 
                    pressure_inc(n) = (2/D*fTP/RHO_f*G^2*(1+...
                        quality_p(i,j,k,m,n)*V_avg*RHO_f) + G^2*V_avg*...
                        (pi*D*q_flux)/(w*H_fv) + G*sin(THETA)/(1/RHO_f+...
                        quality_p(i,j,k,m,n)*V_avg))*(z_increment);
                end
                pressure_drop(i,j,k,m) = sum(pressure_inc);
            end
        end
    end
end

pressure_min = min(pressure_drop);
pressure_max = max(pressure_drop);

