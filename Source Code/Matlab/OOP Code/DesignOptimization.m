% ------------ Design Optimization Routine ------------------
%                    Nicholas Evich
% -----------------------------------------------------------

% Some (or all) of these constraints belong in the actual class definitions

% Minimum and maximum channel diameters (m) as outlined in [Kandlikar, 1990]
d_min = 0.005; % via [Perroud, 1960]
d_max = 0.032; % via [Morzov, 1969]

% Minimum and maximum heat fluxes (W/m^2) as outlined in [Kandlikar, 1990]
qflux_min = 4.7e3; % via [Wright, 1961]
qflux_max = 933e3; % via [Stone, 1971]
% qflux_max = 2280e3; % via [Perroud, 1960]

% Minimum and maximum flow qualities as outlined in [Kandlikar, 1990]
x_min = 0.001; % via [Wright, 1961]
x_max = 0.671; % via [Kenning, 1987]
% x_max = 0.699; % via [Perroud, 1960]
% Maximum flow quality for heat transfer calculations []
% x_max = 0.8; []

% Minimum and maximum mass fluxes (kg/m^s*s) as outlined in [Kandlikar, 1990]
G_min = 67; % via [Stone, 1971]
G_max = 2434; % via [Wright, 1961]
% G_max = 8179; % via [Perroud, 1960]
