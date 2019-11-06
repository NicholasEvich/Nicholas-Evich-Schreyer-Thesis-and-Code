%

clc, clear;

% ---------------- Defining all constants-------------------------- %
% --------- ALL CONSTANTS BEGIN WITH A CAPITAL LETTER ------------- %
% Length of pipe (m)
L = 1;
% Gravitational constant (m/s^2)
Grav = 9.81;
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
% Difference between specific volumes of saturated water/steam at 373.15 K (m^3/kg)
V_diff = 1 / RHO_v - 1 / RHO_f;
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
pipe_angle = linspace(0, 45, 100); 
% heat_input = linspace(0,500,30);
heat_input = 500000;
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
% h_inc = linspace(0,0,250); % For 1-D
h_inc = zeros(length(massflow), length(pipe_angle), length(inlet_diameter), length(heat_input), length(z));
h_inc1 = zeros(length(massflow), length(pipe_angle), length(inlet_diameter), length(heat_input), length(z));
h_inc2 = zeros(length(massflow), length(pipe_angle), length(inlet_diameter), length(heat_input), length(z));
h = zeros(length(massflow), length(pipe_angle), length(inlet_diameter), length(heat_input));
% pressure_inc = linspace(0, 0, 249); % Why is this necessary
pressure_inc = zeros(length(massflow), length(pipe_angle), length(inlet_diameter), length(heat_input), length(z) - 1);
pressure_drop = zeros(length(massflow), length(pipe_angle), length(inlet_diameter), length(heat_input));
hsp = zeros(length(massflow), length(pipe_angle), length(inlet_diameter), length(heat_input), length(z));
FFr = zeros(length(massflow), length(pipe_angle), length(inlet_diameter), length(heat_input), length(z));


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
                       error('Flow quality too high (greater than 0.8). Heat transfer coefficient calculations beyond this point are invalid');
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
                    hsp(i,j,k,m,n) = (K_f/D)*(darcyf/8*(Re - 1000)*PR)/...
                        (1 + 12.7*(darcyf/8)^0.5*(PR^(2/3)-1)); % Calculating
                    % single phase convection coefficient
                    
                    % Ensure that this is supposed to be grav and not G
                    % (areal flow rate)
                    if (w/(A*RHO_f))^2/(Grav*D) >= 0.04 % Conditional statement to 
                        % determine FFr (what is FFr)
                        FFr(i,j,k,m,n) = 1;
                    else % Figure out where this equation comes from and make s
                        % sure that the order of operations is working
                        % correctly
                        FFr(i,j,k,m,n) = 2.63*((w/(A*RHO_f))^2/(Grav*D))^0.3;
                    end
                    
                    % Heat transfer coefficient increment eq 1 (check
                    % OoO's)
                    h_inc1(i,j,k,m,n) = hsp(i,j,k,m,n)*(0.6683*(RHO_f/RHO_v)^0.1*...
                        quality_h(i,j,k,m,n)^0.16*...
                        (1-quality_h(i,j,k,m,n))^0.64...
                        *FFr(i,j,k,m,n)+1058*(q_flux*A/w/H_fv)^0.7*(1-... % Whats up with all these divisions?
                        quality_h(i,j,k,m,n))^0.8*G_sf);
                    % Heat transfer coefficient increment eq 2
                    h_inc2(i,j,k,m,n) = hsp(i,j,k,m,n)*(1.136*(RHO_f/RHO_v)^0.45*...
                        quality_h(i,j,k,m,n)^0.72*...
                        (1-quality_h(i,j,k,m,n))^0.08...
                        *FFr(i,j,k,m,n)+667.2*(q_flux*A/w/H_fv)^0.7*(1-... % Whats up with all these divisions?
                        quality_h(i,j,k,m,n))^0.8*G_sf);
                    if h_inc1(i,j,k,m,n) >= h_inc2(i,j,k,m,n) 
                        h_inc(i,j,k,m,n) = h_inc1(i,j,k,m,n);
                        % h_inc(n) = h_inc1; % For 1-D
                    else
                        h_inc(i,j,k,m,n) = h_inc2(i,j,k,m,n);
                        % h_inc(n) = h_inc2; % for 1-D
                    end
                    
                    % Why is this here?
                    if sqrt(SIGMA/(Grav*(RHO_f-RHO_v)))/D >= 0.5 % Setting
                        % confinement number limit
                        error('Confinement number too high (greater than 0.5). Heat transfer coefficient calculations beyond this point are invalid');
                    else
                        h(i,j,k,m) = mean(h_inc(i,j,k,m,:));
                        % h(i,j,k,m) = mean(h_inc); % for 1-D
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
                    % pressure_inc(i,j,k,m,n)
                    pressure_inc(i,j,k,m,n) = (2/D*fTP/RHO_f*G^2*(1+...
                        quality_p(i,j,k,m,n)*V_diff*RHO_f) + G^2*V_diff*...
                        (pi*D*q_flux)/(w*H_fv) + Grav*sind(THETA)/(1/RHO_f+... % ERROR: G is not g
                        quality_p(i,j,k,m,n)*V_diff))*(z_increment);
                end
                pressure_drop(i,j,k,m) = sum(pressure_inc(i,j,k,m,:));
                % pressure_drop(i,j,k,m) = sum(pressure_inc);
            end
        end
    end
end

h_inverse = 1./h;

% For the 4 calculations below, I may need to add the (:) to each of the 4D
% arrays when I add more inputs
[pressure_min, pressure_min_index] = min(pressure_drop);
[pressure_max, pressure_max_index] = max(pressure_drop);
[h_inv_min, h_inv_min_index] = min(h_inverse);
[h_inv_max, h_inv_max_index] = max(h_inverse);

% Normalized pressure and heat transfer vector
normalized_pressure = (pressure_drop - pressure_min)./(pressure_max - pressure_min);
normalized_h_inv = (h_inverse - h_inv_min)./(h_inv_max - h_inv_min);

distance = sqrt(normalized_pressure.*normalized_pressure + normalized_h_inv.*normalized_h_inv);

% Removing the : from distance(_) below ___
[distance_min, distance_min_index] = min(distance);
[h_max, h_max_index] = max(h); % Potentially unnecessary line

% Enbsure that this is the proper use of these functions
% Yes, all of these indices should ultimately be 4D arrays (such that I can
% get optimizations for every possible input), but I will need to process
% them further when I am not using all of these input parameters
[a,b,c,d] = ind2sub(size(distance), distance_min_index); % Dont necessarily be concerned about getting different indices for this, they should be different
[a_h,b_h,c_h,d_h] = ind2sub(size(h), h_max_index);
[a_p,b_p,c_p,d_p] = ind2sub(size(pressure_drop), pressure_min_index);

% Find a more efficient way to to these (use matrices)
massflow_utopia = massflow(a);
angle_utopia = pipe_angle(b);
diameter_utopia = inlet_diameter(c);
heat_utopia = heat_input(d);

massflow_max_h = massflow(a_h);
angle_max_h = pipe_angle(b_h);
diameter_max_h = inlet_diameter(c_h);
heat_max_h = heat_input(d_h);

massflow_min_pressure = massflow(a_p);
angle_min_pressure = pipe_angle(b_p);
diameter_min_pressure = inlet_diameter(c_p);
heat_min_pressure = heat_input(d_p);

%
% Something is wrong with this variable
% exit_quality = quality(1,reshape(distance_min_index(1,:,1,:),4,[]),1,:,250);

