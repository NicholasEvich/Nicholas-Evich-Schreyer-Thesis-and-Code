clc, clear;
% Input ranges
% massflow = linspace(0.5,1.5,10); % Mass flow rate (kg/s)
massflow = 0.5;
% Inlet_Diameter = linspace(0.006,0.1, 10); % Inlet diameter satisfying confinement number (m)
Inlet_Diameter = 0.01;
Pipe_Angle = linspace(0,45,100); % (degrees)
% heat_input = linspace(5, 5000000, 10); % Heat input (W)
heat_input = 500;

% Constants
L = 1; % Length of pipe
z = linspace(0, L, 250); % Breaking the emtire pipe length into smaller subsections
g = 9.81; % Gravitatioinal constant (m/s^2)
% For the incline below, this is for work done by gravity. not the angle
% associated with the changing pipe diameter
theta = 0; % Incline of pipe (degrees)
Ti = 373.15; % Inlet temperature of saturated water (pure water?)
Tsat = 373.15; % Saturation temperature of water (K)
hfg = 2256e3; % Latent heat of vaporization (J/kg)
rhof = 958.05; % Density of saturated water at 100 deg C(kg/m^3)
rhog = 0.590; % Density of saturated steam at 100 deg C (kg/m^3)
Cp_f = 4181.3; % Specific heat of liquid water at constant pressure (J/kg*K)
k_f = 0.67703; % Thermal conductivity of water at Tsat (W/mK)
mu_f = 0.0002814; % Dynamic viscosity of water at Tsat (Pa-s)
mu_g = 11.97e-6; % Dynamic viscosity of steam at Tsat (Pa-s)
Pr = 1.76; % Prandlt number of saturated water at Tsat
vfg = 1/rhog - 1/rhof; % Average specific volume of water and gas (m^3/kg)
sigma = 58.9e-3; % Surface tension of water at Tsat (N/m)
Gsf = 1; % Factor for calculating h for water

% YOU COULD LITERALLY REPLACE COUNT1,2,3,4,5, ETC WITH THE W, PHI, ETC.
%   OTHERWISE IT IS THE BIGGEST WASTE OF SPACE AND TIME
% Calculating vapor quality, void fraction, heat transfer coefficient, and
% pressure drop

% ALL OF THESE NESTED LOOPS ARE NOT NECESSARY, 3 OF THE FIRST 4 ARE
% "VECTORS" WITH ONE ELEMENT
% This does raise the point that its difficult to maximize 3 or 4 variables
% at once, this is why (maybe surface topology) and Tagucci array analysis
% is helpful and should be explored (later on). Explore 3D plots as well.
for count1 = 1:numel(massflow) % "count1" sounds like a toddler made it. Use i, j, k, x, and y
% Looping through mass flow rates
    W = massflow(count1); % Move initialization out of loop?
    % Giving mass flow rate variable name ^
    for count2 = 1:numel(Pipe_Angle) % Looping through pipe angle
        phi = Pipe_Angle(count2); % Giving pipe angle name "phi"
        for count3 = 1:numel(Inlet_Diameter) % Looping through diameters
            % Both of these are a waste of space. I could literally make
            %   A_sur my counter variable
            % Also, I need to fix this area thing because its not very
            %   useful (see heat flux)
            % Verify this equation for the pipe is correct
            % Figure out if the pipe angles in or out and try both
            % Keep in mind that a contracting area will cause a pressure
            %   drop and a decrease in heat flux (because its a smaller area)
            % Try two things: Keeping the heat flux constant over each
            %   subsection of the pipe and keeping the heat transfer (units
            %   W) constant over each section of the pipe. See which gives
            %   more interesting results
            Do = Inlet_Diameter(count3); % Giving inlet diameter variable name "Do"
            A_sur = pi*Do*L + pi*L^2*tand(phi)/2; % Calculating surface area of pipe
            for count4 = 1:numel(heat_input) % Looping through (total) heat input
                q = heat_input(count4)/A_sur; % Giving heat input variable
                % the name "q". This is the HEAT FLUX, and this is a bad
                % way of calculating it. Make the heat flux constant across
                % the entire pipe. 
                % Run four tests:
                %       1.) Pipe angled in, constant heat input
                %       2.) Pipe angled in, constant heat flux
                %       3.) Pipe angled out, constant heat input
                %       4.) Pipe angled out, constant heat flux
                % Keep in mind, a constant heat flux would make more
                % logical sense for most applications
                for count5 = 1:numel(z) % Breaking the pipe up into z intervals
                    z_incl = z(count5); % Giving each z interval a name
                    if z_incl < L/2 % Calculating diameter as a function of axial distance
                        % I'm pretty sure this is z dependent too, because
                        % diameter literally is axial distance
                        D = Do + 2*z(count5)*tand(phi); % tand is tangent in degrees
                    else
                        D = Do + 2*L*tand(phi) - 2*z(count5)*tand(phi);
                        % Make sure these calculations are actually correct
                    end
                    A = (pi/4)*(D^2); % Calculating area from diameter
                    
                    % Flow quality (these calculations may be incorrect, ow
                    % they are being influenced incorrectly by the heat
                    % flux
                    % Preallocate space for Xe
                    Xe(count1, count2, count3, count4, count5) = ...
                        (-Cp_f*(Tsat-Ti))/hfg + (pi*D*q*z(count5))/(W*hfg);
                    
                    % There is definitely a better way to do the below
                    % statement
                    % Make functional .m files
                    if Xe(count1, count2, count3, count4, count5) > 1 % Cutting off 
                        % quality when there is 100% vapor
                        Xe(count1, count2, count3, count4, count5) = 1;
                    end
                    
                    % Void fraction
                    alpha(count1, count2, count3, count4, count5) = ...
                        1/(1 + rhog/rhof)*...
                        (1-Xe(count1,count2,count3,count4,count5));
                    % Calculating void fraction
                    % If I need to type out Xe(count1, count2... more than
                    % once, Im doing something wrong
                end
                
                
      % The below statement is incorrect. Nick did not recalculate D before
      % each calculation of x for heat transfer, so its resusing D from
      % loop/count 5
                % Heat Transfer Coefficient
                for count6 = 1:numel(z) % Breaking up the pipe
                    % Completely unnecessary memory allocation, just reset
                    % and reuse count5 or just use the vector z
                    
                    % ----------------------------------------------------------
                    z_inc2 = z(count6); % Giving each z increment a name
                    if z_inc2 < L/2 %Calculating diameter as a function of axial
                        % pipe distance (YOU LITERALLY JUST DID THIS DUDE
                        % WHY DO IT AGAIN)
                        % Check this fow when zInc actually equals L/2
                        D = Do + 2*z(count6)*tand(phi);
                    else
                        % Supposedly this number is right, I really dont
                        % know how
                        D = Do + 2*L*tand(phi) - 2*z(count6)*tand(phi);
                    end
                    % -----------------------------------------------------------
                    
                    Xeh(count1, count2, count3, count4, count6) = ...
                        (-Cp_f*(Tsat-Ti))/hfg + (pi*D*q*z(count6))/(W*hfg);
                    % Calculating vapor quality for heat transfer
                    % coefficient
                    if Xeh(count1, count2, count3, count4, count6) > 0.8
                        % Setting limit of heat transfer equation at vapor
                        % quality < 0.8
                        Xeh(count1, count2, count3, count4, count6) = NaN;
                    end
                    
                    % -----------------------------------------------------------
%                     z_inc2 = z(count6); % Giving each z increment a name
%                     if z_inc2 < L/2 %Calculating diameter as a function of axial
%                         % pipe distance (YOU LITERALLY JUST DID THIS DUDE
%                         % WHY DO IT AGAIN)
%                         % Check this fow when zInc actually equals L/2
%                         D = Do + 2*z(count6)*tand(phi);
%                     else
%                         % Supposedly this number is right, I really dont
%                         % know how
%                         D = Do + 2*L*tand(phi) - 2*z(count6)*tand(phi);
%                     end
                    % -------------------------------------------------------------
                    
                    A = (pi/4)*D^2; % Hoist this to global scope, or at least a 
                    % level up so that I dont have to do it again
                    
                    u = W/(rhof.*A); % Calculating velocity of single phase (m/s)
                    Re = rhof*u*D/mu_f; % Reynolds number at specific conditions
                    
                   
                    % Setting values for laminar and too turbulent input
                    % conditions
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
                    hsp = (k_f/D)*(darcyf/8*(Re - 1000)*Pr)/...
                        (1 + 12.7*(darcyf/8)^0.5*(Pr^(2/3)-1)); % Calculating
                    % single phase convection coefficient
                    
                    if (W/(A*rhof))^2/(g*D) >= 0.04 % Conditional statement to 
                        % determine FFr (what is FFr)
                        FFr = 1;
                    else % Figure out where this equation comes from and make s
                        % sure that the order of operations is working
                        % correctly
                        FFr = 2.63*((W/(A*rhof))^2/(g*D))^0.3;
                    end
                    % Heat transfer coefficient increment eq 1 (check
                    % OoO's)
                    h_inc1 = hsp*(0.6683*(rhof/rhog)^0.1*...
                        Xeh(count1,count2,count3,count4,count6)^0.16*...
                        (1-Xeh(count1,count2,count3,count4,count6))^0.64...
                        *FFr+1058*(q*A/W/hfg)^0.7*(1-...
                        Xeh(count1,count2,count3,count4,count6))^0.8*Gsf);
                    % Heat transfer coefficient increment eq 2
                    h_inc2 = hsp*(1.136*(rhof/rhog)^0.45*...
                        Xeh(count1,count2,count3,count4,count6)^0.72*...
                        (1-Xeh(count1,count2,count3,count4,count6))^0.08...
                        *FFr+667.2*(q*A/W/hfg)^0.7*(1-...
                        Xeh(count1,count2,count3,count4,count6))^0.8*Gsf);
                    if h_inc1 >= h_inc2 % Choosing which h equation to use
                        h_inc(count6) = h_inc1;
                    else
                        h_inc(count6) = h_inc2;
                    end
                    
                    if sqrt(sigma/(g*(rhof-rhog)))/D >= 0.5 % Setting
                        % confinement number limit
                        h(count1,count2,count3,count4) = NaN;
                    else
                        h(count1,count2,count3,count4) = mean(h_inc);
                    end
                end
                
                % Pressure drop
                    % I really dont think this loop is supposed to be
                    % inside of loop 5. It would make sense for it to be on
                    % the same level as loop 5 and 6
                    for count7 = 1:(numel(z) - 1)
                        % Calculating vapor quality for pressure drop (why
                        % is this any different than the vapor quality in
                        % the other 10000 loops?)
                        % Also, why is the index for Z count7 + 1?
                        % WHY IS COUNT6 REFERENCED? that is wrong on so
                        % many levels, because at this point in the code,
                        % count6 can be referenced but wont change (it will
                        % be at its maximum value)
                        %       Figure how this affects Xep (negatively)
                        %       and fix it
                        % Note this to Dr. Rau, as this could have been
                        % leading to the exit vapor tquality errors
                        
                        
                        % Again, the error here is that vapor quality is
                        % being calculated before D
                        
                        % -------------------------------------------------------------------------------
                        z_inc3 = z(count7); % Giving each z increment a name
                        if z_inc3 < L/2
                            % Again, we already calculated this, dont do it
                            % again... just initialize all of this outside
                            % of the loop
                            D = Do + 2*z(count7)*tand(phi);
                        else
                            D = Do + 2*L*tand(phi) - 2*z(count7)*tand(phi); % THIS IS THE LINE THAT IS MAKING D A VECTOR
                        end
                        % -------------------------------------------------------------------------------
                        
                        Xep(count1,count2,count3,count4,count7) = (-Cp_f*... %NOTE THAT I SWITCHED COUNT6 TO COUNT7 IN THIS LINE TO MAKE THE CODE WORK
                            (Tsat-Ti))/hfg+(pi*D*q*z(count7 + 1))/(W*hfg);
                        if Xep(count1,count2,count3,count4,count7) > 1
                            % Max vapor quality cutoff
                            Xep(count1,count2,count3,count4,count7) = 1;
                        end
                        
                        % ----------------------------------------------------------------------------------------
%                         z_inc3 = z(count7); % Giving each z increment a name
%                         if z_inc3 < L/2
%                             % Again, we already calculated this, dont do it
%                             % again... just initialize all of this outside
%                             % of the loop
%                             D = Do + 2*z(count7)*tand(phi);
%                         else
%                             D = Do + 2*L*tand(phi) - 2*z(count7)*tand(phi); % THIS IS THE LINE THAT IS MAKING D A VECTOR
%                         end
                        % ----------------------------------------------------------------------------------------
                        
                        A = pi/4*D^2; % Calculating area from diameter
                        G = W/A; % Calculating mass flow rate per area(kg/m^2*s)
                        % Isnt mass flow a vector? How is that supposed to
                        % work?
                        if G*D/mu_f < 2300 % Determining c and n constants for 
                            %two-phase friction factor
                            c = 16;
                            n = 1;
                        elseif G*D/mu_f > 20000
                            c = 0.046;
                            n = 0.20;
                        else
                            c = 0.079;
                            n = 0.25;
                        end
                        ffo = c/(G*D/mu_f)^n; % Determining liquid only friction factor
                        % Two-phase average viscosity
                        mu_bar = Xep(count1,count2,count3,count4,count7)*mu_g...
                            + (1-Xep(count1,count2,count3,count4, count7))*mu_f;
                        % Cicchitti two phase friction factor calculation
                        fTP = ffo*(mu_bar/mu_f)^n;
                        % Pressure drop increments through pipe (Pa)
                        pressure_inc(count7) = (2/D*fTP/rhof*G^2*(1+...
                            Xep(count1,count2,count3,count4,count7)*vfg*...
                            rhof) + G^2*vfg*(pi*D*q)/(W*hfg) + g*sin(theta)...
                            /(1/rhof+Xep(count1,count2,count3,count4,count7)...
                            *vfg))*(z(count7+1)-z(count7));
                    end
                    %Summing pressure drop increments to determine total pressure drop
                    pressure_drop(count1,count2,count3,count4) = sum(pressure_inc);
            end
        end
    end
end

% Normalization of pressure drop and heat transfer inverse
% None of the colons are necessary, rather, they are implied
Pmin = min(pressure_drop(:)); % How are these supposed to work when pressure_drop is actually a 4D array?
% is the colon necessary then?
Pmax = max(pressure_drop(:));
h_inversemin = min(1./h(:));
h_inversemax = max(1./h(:));

NormalizedPressure = (pressure_drop - Pmin)./(Pmax - Pmin);

% Why do it like this lmao
Norm_pressure_1_d = NormalizedPressure(:); % Creating a 1D matrix of normalized pressure drop
% Calculating 1D matric of normalized inverse heat transfer coefficient
% (ok, i guess...)
NormalizedHeatTransferInverse = (1./h - h_inversemin)/(h_inversemax - h_inversemin);

% Displaying maximum heat transfer and minimum pressure drop indices
dist_heat_pressure = sqrt(NormalizedHeatTransferInverse.^2 + NormalizedPressure.^2);
% Finding utopia point indices
[mindist_heat_pressure, minidx] = min(dist_heat_pressure(:)); % Im pretty sure I can leave these :'s out of the expression
[i, j, k, l] = ind2sub(size(dist_heat_pressure), minidx); % What does this function do?
exit_quality = Xe(i, j, k, l, numel(z));
fprintf('The input conditions for the utopia point are mass flow = %3.4f (kg/s), pipe angle = %3.4f (deg),\ninlet diameter = %3.4f (m), and heat input = %3.4f (W).\nThe utopia exit flow quality = %3.4f.\n',massflow(i),Pipe_Angle(j),Inlet_Diameter(k),heat_input(l),exit_quality)

[max_heat, maxidx] = max(h(:)); %Finding highest heat transfer indices
[m, p, r, s] = ind2sub(size(h),maxidx);
exit_quality = Xe(m,p,r,s,numel(z)); % This is so bad, there's no reason to program like this
fprintf('The highest heat transfer has mass flow = %3.4f (kg/s), pipe angle = %3.4f (deg),\ninlet diameter = %3.4f (m), and heat input = %3.4f (W).\nThe exit flow quality = %3.4f.\n',massflow(m),Pipe_Angle(p),Inlet_Diameter(r),heat_input(s),exit_quality)

[min_pressure, minidx] = min(pressure_drop(:)); %Finding lowest pressure drop indices
[t, v, w, y] = ind2sub(size(pressure_drop),minidx);
exit_quality = Xe(t,v,w,y,numel(z));
fprintf('The lowest pressure drop has mass flow = %3.4f (kg/s), pipe angle = %3.4f (deg),\ninlet diameter = %3.4f (m), and heat input = %3.4f (W).\nThe exit flow quality = %3.4f.\n',massflow(t),Pipe_Angle(v),Inlet_Diameter(w),heat_input(y),exit_quality)