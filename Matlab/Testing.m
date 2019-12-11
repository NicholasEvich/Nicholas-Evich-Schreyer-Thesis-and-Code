%% Test for heat transfer coefficient equation

clc;clear

x_low = 0:0.001:0.08;
x_high = 0:0.01:0.8;

term1_low = x_low.^(0.16).*(1-x_low).*(0.64);
term2_low = x_low.^(0.72).*(1-x_low).*(0.08);

term1_high = x_high.^(0.16).*(1-x_high).*(0.64);
term2_high = x_high.^(0.72).*(1-x_high).*(0.08);

figure;
plot(x_low, term1_low);
hold on;
plot(x_low, term2_low);
title('Low Heat Input');
xlabel('Flow Quality');
ylabel('Selected Term from h Equation');
legend('h term 1', 'h term 2');

figure;
plot(x_high, term1_high);
hold on;
plot(x_high, term2_high);
title('High Heat Input');
xlabel('Flow Quality');
ylabel('Selected Term from h Equation');
legend('h term 1', 'h term 2');

%% Test for full version of the pressure drop equation (11/24/19)
clear;clc

W = 0.5;
D = 0.1;
A = 0.25*pi*D*D;
P0 = 100000;
pressure = zeros(1,1000);
pressure_drop = zeros(1,999);
v_f = 1/958.05;
v_g = 1/0.590;
v_fg = v_g - v_f;
h_fg = 2256E03;
q_flux = 1e7;
G = W/A;
MU_f = 0.0002814;
MU_v = 11.97E-06;
theta = 0;

x = linspace(0,0.8,1000);

KE = zeros(1, length(x));
Fric = zeros(1, length(x));
Acc = zeros(1, length(x));
Grav = zeros(1, length(x));
Compress = zeros(1, length(x));
Flash = zeros(1, length(x));

pressure(1) = P0;
L = 1;
dz = L/length(x);

for j = 1:length(x)
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
    MU_bar = x(j)*MU_v...
        + (1-x(j))*MU_f;
    fTP = ffo*(MU_bar/MU_f)^N;
    
    KE(j) = KE_PhaseChange(G, h_fg, v_f, v_g, v_fg, x(j));
    Fric(j) = dP_Friction(D, fTP, G, v_f, v_fg, x(j));
    Grav(j) = dP_Gravity(theta, v_f, v_fg, x(j));
    Acc(j) = dP_Acceleration(D, G, h_fg, q_flux, v_fg, W);
    Compress(j) = 1e-06*Compressibility(G, pressure(j), x(j));
    Flash(j) = 1e-06*Flashing(G, h_fg, pressure(j), v_fg, x(j));
    
    pressure_drop(j) = dz*(KE(j)*Fric(j) + Acc(j) + Grav(j))/(KE(j) + Compress(j) - Flash(j));
    pressure(j + 1) = pressure(j) - pressure_drop(j);
    
end

%% Test for void fraction plot shape (12/2/19)

% RHO_v = 0.590;
% RHO_f = 958.05;
P = linspace(0.101330, 0.050000, 1000);
z = linspace(0,1,1000);
x = linspace(0,1,1000);
v_g = 0.1909.*(P).^(-0.96);
rho_g = 1./v_g;
v_f = - 3e-07.*(P).^6 + 4e-06.*(P).^5 - 3e-05.*(P).^4 + 1e-04.*(P).^3 - 0.0002.*(P).^2 + 0.0002.*(P) + 0.001;
rho_f = 1./v_f;
alpha = 1./(1 + rho_g./rho_f.*(1 - x)./x);


plot(z, alpha);