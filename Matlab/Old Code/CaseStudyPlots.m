%% -----------------Straight Pipe, No Heat Input--------------------------

% Flow Quality
%   There should be no real trend here
plot(z, reshape(quality(1,1,1,1,:), 1, []));
xlabel('Axial Distance (m)');
ylabel('Instantaneous Flow Quality');

% Pressure Drop
%   Incremental Pressure Drop
figure;
plot(z(2:length(z)), reshape(pressure_inc(1,1,1,1,:), 1, []));
xlabel('Axial Distance (m)');
ylabel('Incremental Pressure Drop (Pa)');
%   Cumulative Pressure Drop
figure;
plot(z(2:length(z)), cumsum(reshape(pressure_inc(1,1,1,1,:), 1, [])));
xlabel('Axial Distance (m)');
ylabel('Cumulative Pressure Drop (Pa)');

% Heat Transfer Coefficient
%   There should be no real trend here
figure;
plot(z, reshape(h_inc(1,1,1,1,:), 1, []));
xlabel('Axial Distance (m)');
ylabel('Incremental Heat Transfer Coefficient');

%% --------------Straight Pipe, Variable Heat Input-----------------------

%% Low Heat Input

% Flow Quality
%   There should be no real trend here
plot(z, reshape(quality(1,1,1,1,:), 1, []));
xlabel('Axial Distance (m)');
ylabel('Instantaneous Flow Quality');

% Pressure Drop
%   Incremental Pressure Drop
figure;
plot(z(2:length(z)), reshape(pressure_inc(1,1,1,1,:), 1, []));
xlabel('Axial Distance (m)');
ylabel('Incremental Pressure Drop (Pa)');
%   Cumulative Pressure Drop
figure;
plot(z(2:length(z)), cumsum(reshape(pressure_inc(1,1,1,1,:), 1, [])));
xlabel('Axial Distance (m)');
ylabel('Cumulative Pressure Drop (Pa)');

% Heat Transfer Coefficient
%   There should be no real trend here
figure;
plot(z, reshape(h_inc(1,1,1,1,:), 1, []));
xlabel('Axial Distance (m)');
ylabel('Incremental Heat Transfer Coefficient');

% Void Fraction
figure;
plot(z, reshape(alpha(1,1,1,1,:), 1, []));
xlabel('Axial Distance (m)');
ylabel('Void Fraction');

%% High Heat Input

% Flow Quality
%   There should be no real trend here
plot(z, reshape(quality(1,1,1,1,:), 1, []));
xlabel('Axial Distance (m)');
ylabel('Instantaneous Flow Quality');

% Pressure Drop
%   Incremental Pressure Drop
figure;
plot(z(2:length(z)), reshape(pressure_inc(1,1,1,1,:), 1, []));
xlabel('Axial Distance (m)');
ylabel('Incremental Pressure Drop (Pa)');
%   Cumulative Pressure Drop
figure;
plot(z(2:length(z)), cumsum(reshape(pressure_inc(1,1,1,1,:), 1, [])));
xlabel('Axial Distance (m)');
ylabel('Cumulative Pressure Drop (Pa)');

% Heat Transfer Coefficient
%   There should be no real trend here
figure;
plot(z, reshape(h_inc(1,1,1,1,:), 1, []));
xlabel('Axial Distance (m)');
ylabel('Incremental Heat Transfer Coefficient');

% Void Fraction
figure;
plot(z, reshape(alpha(1,1,1,1,:), 1, []));
hold on;
plot(z, reshape(quality(1,1,1,1,:), 1, []));
xlabel('Axial Distance (m)');
ylabel('Void Fraction');

%% Trends Across Heat Input

% Flow Quality
%   Flow Quality vs Axial Distance
plot(z, reshape(quality(1,1,1,1,:), 1, []));
hold on;
for numHeat = 2:length(heat_input)
    plot(z, reshape(quality(1,1,1,numHeat,:), 1, []));
end
xlabel('Axial Distance (m)');
ylabel('Instantaneous Flow Quality');
%   Exit Flow Quality for Many Heat Inputs
figure;
plot(heat_input, reshape(quality(1,1,1,:,length(z)), 1, []));
xlabel('Heat Input (W)');
ylabel('Exit Flow Quality');

% Pressure Drop
%   Incremental Pressure Drop
figure;
plot(z(2:length(z)), reshape(pressure_inc(1,1,1,1,:), 1, []));
hold on;
for numHeat = 2:length(heat_input)
    plot(z(2:length(z)), reshape(pressure_inc(1,1,1,numHeat,:), 1, []));
end
xlabel('Axial Distance (m)');
ylabel('Incremental Pressure Drop (Pa)');
%   Cumulative Pressure Drop
figure;
plot(z(2:length(z)), cumsum(reshape(pressure_inc(1,1,1,1,:), 1, [])));
hold on;
for numHeat = 2:length(heat_input)
    plot(z(2:length(z)), cumsum(reshape(pressure_inc(1,1,1,numHeat,:), 1, [])));
end
xlabel('Axial Distance (m)');
ylabel('Cumulative Pressure Drop (Pa)');
%   Total Pressure Drop for Many Heat Inputs
figure;
plot(heat_input, reshape(pressure_drop(1,1,1,:), 1, []));
xlabel('Heat Input (W)');
ylabel('Total Pressure Drop (Pa)');

% Heat Transfer Coefficient
%   Incremental Heat Transfer Coefficient
%       TBD
%   Overall Heat Transfer Coefficients for Many Heat Inputs
figure;
plot(heat_input, reshape(h(1,1,1,:), 1, []));
xlabel('Heat Input (W)');
ylabel('Overall Heat Transfer Coefficient');
%% ------------------Angled Pipe, No Heat Input---------------------------

% Flow Quality
%   There should be no real trend here
plot(z, reshape(quality(1,5,1,1,:), 1, []));
xlabel('Axial Distance (m)');
ylabel('Instantaneous Flow Quality');

% Pressure Drop
%   Incremental Pressure Drop
figure;
semilogy(z(2:length(z)), reshape(pressure_inc(1,1,1,1,:), 1, []));
hold on;
for numAngle = 2:length(pipe_angle)
    semilogy(z(2:length(z)), reshape(pressure_inc(1,numAngle,1,1,:), 1, []));
end
xlabel('Axial Distance (m)');
ylabel('Incremental Pressure Drop (Pa)');
%   Cumulative Pressure Drop
figure;
semilogy(z(2:length(z)), cumsum(reshape(pressure_inc(1,1,1,1,:), 1, [])));
hold on;
for numAngle = 2:length(pipe_angle)
    semilogy(z(2:length(z)), cumsum(reshape(pressure_inc(1,numAngle,1,1,:), 1, [])));
end
xlabel('Axial Distance (m)');
ylabel('Cumulative Pressure Drop (Pa)');

% Heat Transfer Coefficient
%   There should be no real trend here
figure;
plot(z, reshape(h_inc(1,5,1,1,:), 1, []));
xlabel('Axial Distance (m)');
ylabel('Incremental Heat Transfer Coefficient');

%% -----------------Angled Pipe, Variable Heat Input----------------------

%% Low Heat Input

% Flow Quality
%   Flow Quality vs Axial Distance
% Pressure Drop
%   Incremental Pressure Drop
%   Cumulative Pressure Drop
%   Total Pressure Drop
% Heat Transfer Coefficient (IDK)
%   Incremental Heat Transfer Coefficient
%% High Heat Input

% Flow Quality
%   Flow Quality vs Axial Distance
% Pressure Drop
%   Incremental Pressure Drop
%   Cumulative Pressure Drop
%   Total Pressure Drop
% Heat Transfer Coefficient (IDK)
%   Incremental Heat Transfer Coefficient
%% Trends Across Heat Input

% Flow Quality
%   Exit Quality for the Same Heat Input, Different Angles
%   Exit Quality for the Same Angle, Different Heat Inputs
%   3D Plot of Exit Quality for Every Heat Input and Angle
% Pressure Drop
%   Total Pressure Drop for the Same Heat Input, Different Angles
%   Total Pressure Drop for the Same Angle, Different Heat Inputs
%   3D Plot of Total Pressure Drop for Every Heat Input and Angle
% Heat Transfer Coefficient
%   Overall Heat Transfer Coefficient for the Same Heat Input, Different Angles
%   Overall Heat Transfer Coefficient for the Same Angle, Different Heat Inputs
%   3D Plot of Overall Heat Transfer Coefficient for Every Heat Input and Angle