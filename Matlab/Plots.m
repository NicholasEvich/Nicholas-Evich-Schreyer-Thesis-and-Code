%%
% For heat_input = linspace(100, 1000000, 30);
plot(heat_input, angle_utopia);
xlabel('Heat Input (W)');
ylabel('Utopia Angle (degrees)');
title('10/28/19: 100 angles constrained between 2 and 3 degrees');
%% 
% For heat_input = linspace(902300, 902400, 100);
plot(heat_input, angle_utopia);
xticks(linspace(902300,902400,6));
xticklabels({'902.30', '902.32', '902.34', '902.36', '902.38', '902.40'});
xlabel('Heat Input (kW)');
ylabel('Utopia Angle (degrees)');

%% Graphs
% Plots that I want to do
% Quality as a function of axial distance
% Exit quality as a function of pipe angle
% Change in quality for each pipe increment (for a given angle)

% Figure out what this is doing and find a cleaner way to do it
% Plotting quality as a function of axial distance for each optimization
% case
figure;
plot(z, reshape(quality(1,b,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
hold on;
plot(z, reshape(quality(1,b_h,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
plot(z, reshape(quality(1,b_p,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
str1 = sprintf(['Utopia Point (%.2f' char(176) ')'], angle_utopia);
str2 = sprintf(['Maximum h (%.2f' char(176) ')'], angle_max_h);
str3 = sprintf(['Minimum Delta P (%.2f' char(176) ')'], angle_min_pressure);
legend(str1, str2, str3, 'Location', 'northwest');
xlabel('Axial Distance Along Pipe (m)');
ylabel('Flow Quality');

% Plotting total pressure drop as a function of pipe angle (linear scale)
figure;
yyaxis left;
plot(pipe_angle, pressure_drop);
ylabel('Total Pressure Drop (Pa)');
yyaxis right;
semilogy(pipe_angle, pressure_drop);
xlabel('Pipe Expansion Angle (degrees)');
legend('Linear','Logarithmic');
title('10/21/19: 902.37 kW, 0.01 m, 0.5 kg/s');

% Plotting incremental pressure drop as a function of pipe angle
figure;
plot(z(2:250), reshape(pressure_inc(1,b,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
hold on;
plot(z(2:250), reshape(pressure_inc(1,b_h,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
plot(z(2:250), reshape(pressure_inc(1,b_p,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
str1 = sprintf(['Utopia Point (%.2f' char(176) ')'], angle_utopia);
str2 = sprintf(['Maximum h (%.2f' char(176) ')'], angle_max_h);
str3 = sprintf(['Minimum Delta P (%.2f' char(176) ')'], angle_min_pressure);
legend(str1, str2, str3, 'location', 'northwest');
xlabel('Axial Distance Along Pipe (m)');
ylabel('Incremental Pressure Drop (Pa)');

% Logarithmic incremental pressure drop
figure;
semilogy(z(2:250), reshape(pressure_inc(1,b,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
hold on;
semilogy(z(2:250), reshape(pressure_inc(1,b_h,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
semilogy(z(2:250), reshape(pressure_inc(1,b_p,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
str1 = sprintf(['Utopia Point (%.2f' char(176) ')'], angle_utopia);
str2 = sprintf(['Maximum h (%.2f' char(176) ')'], angle_max_h);
str3 = sprintf(['Minimum Delta P (%.2f' char(176) ')'], angle_min_pressure);
legend(str1, str2, str3, 'location', 'northwest');
% title('');
xlabel('Axial Distance Along Pipe (m)');
ylabel('Incremental Pressure Drop (Pa)');

% Linear Cumulative Pressure Drop
figure;
plot(z(2:250), cumsum(reshape(pressure_inc(1,b,1,1,:), 1, []))); % massflow, utopia angle, inlet diam, heat, every z increment
hold on;
plot(z(2:250), cumsum(reshape(pressure_inc(1,b_h,1,1,:), 1, []))); % massflow, utopia angle, inlet diam, heat, every z increment
plot(z(2:250), cumsum(reshape(pressure_inc(1,b_p,1,1,:), 1, []))); % massflow, utopia angle, inlet diam, heat, every z increment
str1 = sprintf(['Utopia Point (%.2f' char(176) ')'], angle_utopia);
str2 = sprintf(['Maximum h (%.2f' char(176) ')'], angle_max_h);
str3 = sprintf(['Minimum Delta P (%.2f' char(176) ')'], angle_min_pressure);
legend(str1, str2, str3, 'location', 'northwest');
% title('');
xlabel('Axial Distance Along Pipe (m)');
ylabel('Cumulative Pressure Drop (Pa)');

% Logarithmic Cumulative Pressure Drop

figure;
semilogy(z(2:250), cumsum(reshape(pressure_inc(1,b,1,1,:), 1, []))); % massflow, utopia angle, inlet diam, heat, every z increment
hold on;
semilogy(z(2:250), cumsum(reshape(pressure_inc(1,b_h,1,1,:), 1, []))); % massflow, utopia angle, inlet diam, heat, every z increment
semilogy(z(2:250), cumsum(reshape(pressure_inc(1,b_p,1,1,:), 1, []))); % massflow, utopia angle, inlet diam, heat, every z increment
str1 = sprintf(['Utopia Point (%.2f' char(176) ')'], angle_utopia);
str2 = sprintf(['Maximum h (%.2f' char(176) ')'], angle_max_h);
str3 = sprintf(['Minimum Delta P (%.2f' char(176) ')'], angle_min_pressure);
legend(str1, str2, str3, 'location', 'northwest');
% title('');
xlabel('Axial Distance Along Pipe (m)');
ylabel('Cumulative Pressure Drop (Pa)');

% Plotting heat transfer coefficient as a function of pipe angle

figure;
yyaxis left;
plot(pipe_angle, h);
ylabel('Heat Transfer Coeffient');
yyaxis right;
semilogy(pipe_angle, h);
xlabel('Pipe Expansion Angle (degrees)');
legend('Linear','Logarithmic');
title('10/21/19: 902.37 kW, 0.01 m, 0.5 kg/s');

% Plotting inverse heat transfer coeffient
figure;
yyaxis left;
plot(pipe_angle, h_inverse);
ylabel('1/h');
yyaxis right;
semilogy(pipe_angle, h_inverse);
xlabel('Pipe Expansion Angle (degrees)');
legend('Linear','Logarithmic');
title('10/21/19: 1 MW, 0.01 m, 0.5 kg/s');

%% Plots for heat transfer coefficient increment

figure;
plot(z, reshape(h_inc(1,b,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
hold on;
plot(z, reshape(h_inc(1,b_h,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
plot(z, reshape(h_inc(1,b_p,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
str1 = sprintf(['Utopia Point (%.2f' char(176) ')'], angle_utopia);
str2 = sprintf(['Maximum h (%.2f' char(176) ')'], angle_max_h);
str3 = sprintf(['Minimum Delta P (%.2f' char(176) ')'], angle_min_pressure);
legend(str1, str2, str3, 'location', 'northwest');
xlabel('Axial Distance Along Pipe (m)');
ylabel('Incremental Heat Transfer Coefficient');
title('h inc');

figure;
semilogy(z, reshape(h_inc(1,b,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
hold on;
semilogy(z, reshape(h_inc(1,b_h,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
semilogy(z, reshape(h_inc(1,b_p,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
str1 = sprintf(['Utopia Point (%.2f' char(176) ')'], angle_utopia);
str2 = sprintf(['Maximum h (%.2f' char(176) ')'], angle_max_h);
str3 = sprintf(['Minimum Delta P (%.2f' char(176) ')'], angle_min_pressure);
legend(str1, str2, str3, 'location', 'northwest');
xlabel('Axial Distance Along Pipe (m)');
ylabel('Incremental Heat Transfer Coefficient');
title('h inc');

figure;
plot(z, reshape(h_inc1(1,b,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
hold on;
plot(z, reshape(h_inc1(1,b_h,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
plot(z, reshape(h_inc1(1,b_p,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
str1 = sprintf(['Utopia Point (%.2f' char(176) ')'], angle_utopia);
str2 = sprintf(['Maximum h (%.2f' char(176) ')'], angle_max_h);
str3 = sprintf(['Minimum Delta P (%.2f' char(176) ')'], angle_min_pressure);
legend(str1, str2, str3, 'location', 'northwest');
xlabel('Axial Distance Along Pipe (m)');
ylabel('Incremental Heat Transfer Coefficient');
title('h inc1');

figure;
plot(z, reshape(h_inc2(1,b,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
hold on;
plot(z, reshape(h_inc2(1,b_h,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
plot(z, reshape(h_inc2(1,b_p,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
str1 = sprintf(['Utopia Point (%.2f' char(176) ')'], angle_utopia);
str2 = sprintf(['Maximum h (%.2f' char(176) ')'], angle_max_h);
str3 = sprintf(['Minimum Delta P (%.2f' char(176) ')'], angle_min_pressure);
legend(str1, str2, str3, 'location', 'northwest');
xlabel('Axial Distance Along Pipe (m)');
ylabel('Incremental Heat Transfer Coefficient');
title('h inc2');

%% Plotting hsp and Froude number (FFr) vs axial distance for the optimization angles
% to determine if they are affecting the value of the heat transfer
% coefficient significantly

figure;
semilogy(z, reshape(hsp(1,b,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
hold on;
semilogy(z, reshape(hsp(1,b_h,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
semilogy(z, reshape(hsp(1,b_p,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
str1 = sprintf(['Utopia Point (%.2f' char(176) ')'], angle_utopia);
str2 = sprintf(['Maximum h (%.2f' char(176) ')'], angle_max_h);
str3 = sprintf(['Minimum Delta P (%.2f' char(176) ')'], angle_min_pressure);
legend(str1, str2, str3, 'location', 'northwest');
xlabel('Axial Distance Along Pipe (m)');
ylabel('hsp');

figure;
semilogy(z, reshape(FFr(1,b,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
hold on;
semilogy(z, reshape(FFr(1,b_h,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
semilogy(z, reshape(FFr(1,b_p,1,1,:), 1, [])); % massflow, utopia angle, inlet diam, heat, every z increment
str1 = sprintf(['Utopia Point (%.2f' char(176) ')'], angle_utopia);
str2 = sprintf(['Maximum h (%.2f' char(176) ')'], angle_max_h);
str3 = sprintf(['Minimum Delta P (%.2f' char(176) ')'], angle_min_pressure);
legend(str1, str2, str3, 'location', 'northwest');
xlabel('Axial Distance Along Pipe (m)');
ylabel('Froude Numbeer (FFr)');

%% Plotting heat transfer coefficient and pressure drop vs pipe angle on the same plots
[h_min, h_min_index] = min(h);
[h_max, h_max_index] = max(h);
normalized_h = (h - h_min)./(h_max - h_min);

figure;
plot(pipe_angle, normalized_h);
hold on;
plot(pipe_angle, normalized_h_inv);
plot(pipe_angle, normalized_pressure);
xlabel('Pipe Expansion Angle');
ylabel('Normalized Output Parameter');
title('10/24/19: 500 W, 0.1 m, 0.5 kg/s');
legend('Heat Transfer Coefficient','Inverse Heat Transfer Coefficient','Pressure Drop');