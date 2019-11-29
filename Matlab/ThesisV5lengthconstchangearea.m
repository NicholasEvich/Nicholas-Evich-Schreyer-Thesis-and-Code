clc, clear
%Input Ranges
%massflow = linspace(0.07,1.5,10); %mass flow rate (kg/s) 
%Inlet_Diameter = linspace(0.006,.1,100); %inlet diameter satisfying
%confinement number (m)
Pipe_Angle = linspace(0,45,100); %pipe expansion angle (degrees)
%heat_input = linspace(5,5000,100); %heat input (W)

%Inputs held constant for two-axes plots
massflow = .5;
Inlet_Diameter = .01;
%Pipe_Angle = 0;
heat_input = 500; %low heat
%heat_input = 35000; %high heat
%heat_input = 0;  %zero heat input

%Constants
L = 1; %length of pipe
z = linspace(0,L,250); %breaking the entire length into smaller subsections
g = 9.81; %gravatational constant (m/s^2)
theta = 0; %incline of pipe (degrees)
Ti = 373.15; %inlet temperature of saturated water (K)
Tsat = 373.15; %water saturation temperature (K)
hfg = 2256*10^3; %latent heat of vaporization (J/kg)
rhof = 958.05; %density of saturated water at 100 deg C (kg/m^3)
rhog = 0.590; %density of saturated steam at 100 deg C (kg/m^3)
Cp_f = 4181.3; %specific heat of liquid water at constant pressure (J/kg*K)
k_f = 0.67703; %thermal conductivity of water at Tsat (W/mK)
mu_f = 0.0002814; %dynamic viscosity of water at Tsat (Pa-s)
mu_g = 11.97E-6; %dynamic viscosity of water vapor at Tsat (Pa-s)
Pr = 1.76; %Prandlt number of saturated water at Tsat
vfg = 1/rhog - 1/rhof; %average specific volume of water and gas (m^3/kg)
sigma = 58.9E-3; %surface tension of water at Tsat (N/m)
Gsf = 1; %factor for calculating h for water

%Calculating vapor quality, void fraction, heat transfer coefficient,
%and pressure drop
for count1 = 1:numel(massflow) %Looping through mass flow rates
    W = massflow(count1); %Giving massflow rate variable name "W"
    for count2 = 1:numel(Pipe_Angle) %Looping through pipe angle
        phi = Pipe_Angle(count2); %Giving pipe angle variable name "phi"
        for count3 = 1:numel(Inlet_Diameter) %Looping through diameters
            Do = Inlet_Diameter(count3); %Giving inlet diameter variable
            %name "Do"
            A_sur = pi*Do*L + pi*L^2*tand(phi)/2; %Calculating suface area
            %of pipe   
            for count4 = 1:numel(heat_input) %Looping through heat input
                q = heat_input(count4)/A_sur; %Giving heat input variable
                %name "q"
                for count5 = 1:numel(z) %Breaking up the pipe into z
                    %intervals
                    z_inc1 = z(count5); %Giving each z interval a name
                    if z_inc1 < L/2 %Calculating diameter as a function of
                        %axial distance
                       D = Do + 2*z(count5)*tand(phi);
                    else
                       D = Do + 2*L*tand(phi) - 2*z(count5)*tand(phi);
                    end
                    A = pi/4*D^2; %Calculating area from diameter
               
                    %Flow quality
                    Xe(count1,count2,count3,count4,count5) = ...
                        (-Cp_f*(Tsat-Ti))/hfg+(pi*D*q*z(count5))/(W*hfg); 
                    %Calculating vapor quality
                    if Xe(count1,count2,count3,count4,count5) > 1 %Cutting
                        %off vapor quality when there is 100% vapor
                       Xe(count1,count2,count3,count4,count5) = 1;
                    end
                    
                    %Void fraction
                    alpha(count1,count2,count3,count4,count5) = ...
                        1/(1+(rhog/rhof)*...
                        (1-Xe(count1,count2,count3,count4,count5))/...
                        Xe(count1,count2,count3,count4,count5)); 
                    %Calculating void fraction
                end
                %Heat transfer coefficient
                for count6 = 1:numel(z) %Breaking up pipe into z intervals 
                    %for heat transfer coefficient calculation
                    Xeh(count1,count2,count3,count4,count6) = ...
                        (-Cp_f*(Tsat-Ti))/hfg+(pi*D*q*z(count6))/(W*hfg); 
                    %Calculating vapor quality for heat transfer
                    %coefficient
                    if Xeh(count1,count2,count3,count4,count6) > 0.8  
                        %Setting limit of heat transfer equations at vapor
                        %quality < 0.8
                       Xeh(count1,count2,count3,count4,count6) = nan;
                    end
                    z_inc2 = z(count6); %Giving each z increment a name
                    if z_inc2 < L/2 %Calculating diameter as a function of
                        %axial pipe distance
                       D = Do + 2*z(count6)*tand(phi);
                    else
                       D = Do + 2*L*tand(phi) - 2*z(count6)*tand(phi);
                    end
                    A = pi/4*D^2; %Calculating area from diameter
                    u = W/(rhof.*A); %Calculating velocity of 
                    %single-phase (m/s)
                    Re = rhof*u*D/mu_f; %Reynolds number at specific 
                    %conditions
                    if Re < 2000 %Setting errors for laminar and too 
                        %turbulent input conditions
                        error('Flow cannot be laminar. Change input conditions, i.e. increase massflow or decrease diameter.')
                    elseif Re > 5E6
                        error('Flow is too turbulent. Change input conditions, i.e. decrease massflow or increase diameter.')
                    end
                    darcyf = (0.79*log(Re)-1.64)^(-2); %Calculating darcy friction factor for smooth pipes
                    hsp = (k_f/D)*(darcyf/8*(Re-1000)*Pr)/(1+12.7*(darcyf/8)^0.5*(Pr^(2/3)-1)); %Calculating single-phase convection coefficient
                    if (W/(A*rhof))^2/(g*D) >= 0.04 %Conditional statement to determine FFr
                        FFr = 1;
                    else
                        FFr = 2.63*((W/(A*rhof))^2/(g*D))^0.3;
                    end
                    h_inc1 = hsp*(0.6683*(rhof/rhog)^0.1*Xeh(count1,count2,count3,count4,count6)^0.16*(1-Xeh(count1,count2,count3,count4,count6))^0.64*FFr+1058*(q*A/W/hfg)^0.7*(1-Xeh(count1,count2,count3,count4,count6))^0.8*Gsf); %Heat transfer coefficient increment eq 1
                    h_inc2 = hsp*(1.136*(rhof/rhog)^0.45*Xeh(count1,count2,count3,count4,count6)^0.72*(1-Xeh(count1,count2,count3,count4,count6))^0.08*FFr+667.2*(q*A/W/hfg)^0.7*(1-Xeh(count1,count2,count3,count4,count6))^0.8*Gsf); %Heat transfer coefficient increment eq 2
                    if h_inc1 >= h_inc2  %Choosing which h equation to use
                       h_inc(count6) = h_inc1;
                    else
                       h_inc(count6) = h_inc2;
                    end
                end
                if sqrt(sigma/(g*(rhof-rhog)))/D >= 0.5 %Setting confinement number limit
                   h(count1,count2,count3,count4) = NaN;
                else 
                   h(count1,count2,count3,count4) = mean(h_inc);
                end
                
                %Pressure drop 
                for count7 = 1:(numel(z)-1) %Breaking up pipe into z - 1 increments for pressure drop calculation
                    Xep(count1,count2,count3,count4,count6) = (-Cp_f*(Tsat-Ti))/hfg+(pi*D*q*z(count7+1))/(W*hfg); %Calculating vapor quality for pressure drop
                    if Xep(count1,count2,count3,count4,count7) > 1   %Condition cutting off vapor quality when there is 100% vapor
                       Xep(count1,count2,count3,count4,count7) = 1;
                    end
                    z_inc3 = z(count7); %Giving each z increment a name
                    if z_inc3 < L/2 %Calculating diameter as a function of axial distance
                       D = Do + 2*z(count7)*tand(phi);
                    else
                       D = Do + 2*L*tand(phi) - 2*z(count7)*tand(phi);
                    end
                    A = pi/4*D^2; %Calculating area from diameter
                    G = W/A; %Calculating mass flow rate per area(kg/m^2*s)
                    if G*D/mu_f < 2300 %Determining c and n constants for two-phase friction factor
                        c = 16;
                        n = 1;
                    elseif G*D/mu_f > 20000
                        c = 0.046;
                        n = 0.20;
                    else 
                        c = 0.079;
                        n = 0.25;
                    end
                    ffo = c/(G*D/mu_f)^n; %Determining liquid only friction factor
                    mu_bar = Xep(count1,count2,count3,count4,count7)*mu_g + (1-Xep(count1,count2,count3,count4,count7))*mu_f; %Calculating two-phase average viscosity
                    fTP = ffo*(mu_bar/mu_f)^n; %Cicchitti two phase friction factor calculation
                    pressure_inc(count7) = (2/D*fTP/rhof*G^2*(1+Xep(count1,count2,count3,count4,count7)*vfg*rhof) + G^2*vfg*(pi*D*q)/(W*hfg) + g*sin(theta)/(1/rhof+Xep(count1,count2,count3,count4,count7)*vfg))*(z(count7+1)-z(count7)); %Calculating pressure drop increment through pipe (Pa)
                end
                pressure_drop(count1,count2,count3,count4) = sum(pressure_inc); %Summing pressure drop increments to determine total pressure drop
            end    
        end    
    end
end

%Normalization of pressure drop and heat transfer inverse
Pmin = min(pressure_drop(:)); %Defining minimum pressure for input ranges
Pmax = max(pressure_drop(:)); %Defining maximum pressure for input ranges
h_inversemin = min(1./h(:)); %Defining minimum for inverse of heat transfer coefficient for input ranges
h_inversemax = max(1./h(:)); %Defining maximum for inverse of heat transfer coefficient for input ranges

NormalizedPressure = (pressure_drop-Pmin)./(Pmax-Pmin); %Calculating normalized pressure drop
Norm_pressure_1_d = NormalizedPressure(:); %Creating 1D matrix of normalized pressure drop
NormalizedHeatTransferInverse = (1./h-h_inversemin)/(h_inversemax-h_inversemin); %Calculating normalized inverse heat transfer coefficient
Norm_h_1_d_inverse = NormalizedHeatTransferInverse(:); %Creating 1D matrix of normalized inverse heat transfer coefficient

%Displaying maximum heat transfer and minimum pressure drop indices
dist_heat_pressure = sqrt(NormalizedHeatTransferInverse.^2 + NormalizedPressure.^2); %Finding utopia point indices
[mindist_heat_pressure, minidx] = min(dist_heat_pressure(:));
[i, j, k, l] = ind2sub(size(dist_heat_pressure), minidx);
exit_quality = Xe(i,j,k,l,numel(z));
fprintf('The input conditions for the utopia point are mass flow = %3.4f (kg/s), pipe angle = %3.4f (deg),\ninlet diameter = %3.4f (m), and heat input = %3.4f (W).\nThe utopia exit flow quality = %3.4f.\n',massflow(i),Pipe_Angle(j),Inlet_Diameter(k),heat_input(l),exit_quality)

[max_heat, maxidx] = max(h(:)); %Finding highest heat transfer indices
[m, p, r, s] = ind2sub(size(h),maxidx);
exit_quality = Xe(m,p,r,s,numel(z));
fprintf('The highest heat transfer has mass flow = %3.4f (kg/s), pipe angle = %3.4f (deg),\ninlet diameter = %3.4f (m), and heat input = %3.4f (W).\nThe exit flow quality = %3.4f.\n',massflow(m),Pipe_Angle(p),Inlet_Diameter(r),heat_input(s),exit_quality)

[min_pressure, minidx] = min(pressure_drop(:)); %Finding lowest pressure drop indices
[t, v, w, y] = ind2sub(size(pressure_drop),minidx);
exit_quality = Xe(t,v,w,y,numel(z));
fprintf('The lowest pressure drop has mass flow = %3.4f (kg/s), pipe angle = %3.4f (deg),\ninlet diameter = %3.4f (m), and heat input = %3.4f (W).\nThe exit flow quality = %3.4f.\n',massflow(t),Pipe_Angle(v),Inlet_Diameter(w),heat_input(y),exit_quality)

%{
%Plot Pressure Drop vs Heat Transfer Coefficient 2 y axis plots
scatter(Norm_h_1_d_inverse,Norm_pressure_1_d)
ylabel('Normalized \DeltaP');
xlabel('Normalized 1/h')
hold on
scatter(NormalizedHeatTransferInverse(i,j,k,l),NormalizedPressure(i,j,k,l),'o','MarkerFaceColor','r')
scatter(NormalizedHeatTransferInverse(m,p,r,s),NormalizedPressure(m,p,r,s),'o','MarkerFaceColor','b')
scatter(NormalizedHeatTransferInverse(t,v,w,y),NormalizedPressure(t,v,w,y),'o','MarkerFaceColor','g')
hold off

xlim([0 .001])
ylim([0 .001])
%}


%{
%Plot Xe vs z
%Xe_1_d = Xe(:);
%alpha_1_d = alpha(:);
Xe_1_d2 = Xe(:);
alpha_1_d2 = alpha(:);
yyaxis right
plot(z,alpha_1_d)
ylabel('Void Fraction')
yyaxis left
plot(z,Xe_1_d)
hold on
yyaxis right
plot(z,alpha_1_d2)
yyaxis left
plot(z,Xe_1_d2)
ylabel('Vapor Quality')
hold off
xlabel('Axial Pipe Distance (m)')
ylim([0 1])
legend('Flow Quality All Vaporized','Flow Quality Not All Vaporized','Void Fraction All Vaporized','Void Fraction Not All Vaporized','location','best')
%}

%{
%Plot 2-Axes Graphs
h_1_d = h(:);
pressure_drop_1_d = pressure_drop(:);
figure
yyaxis left
plot(Pipe_Angle,pressure_drop_1_d)
ylabel('\DeltaP (Pa)')
yyaxis right
plot(Pipe_Angle,h_1_d)
ylabel('h (W/m^2K)')
xlabel('\Phi (\circ)')
%}