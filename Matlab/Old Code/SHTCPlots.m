%% Pareto Plot

norm_p_1d = reshape(normalized_pressure, 1, []);
norm_h_1d = reshape(normalized_h_inv, 1, []);
scatter(norm_p_1d, norm_h_1d);
hold on;
scatter(normalized_pressure(a,b,c,d), normalized_h_inv(a,b,c,d), 500, 'r.');
hold on;
scatter(1,0,500, 'r.');
hold on;
scatter(0,1, 500, 'r.');
xlabel('\DeltaP');
ylabel('1/h');
% set(gca,'FontSize',12);

%% diameter and heat input

heat_input = linspace(500,50000,31);

[distance_min2, distance_min_index2] = min(distance(:,:,:,:));
distance_min2 = reshape(distance_min2, length(inlet_diameter), length(heat_input));
distance_min_index2 = reshape(distance_min_index2, length(inlet_diameter), length(heat_input));

figure;
    plot(0.001*heat_input, pipe_angle(distance_min_index2(13,:)),'b-', 'MarkerSize', 50);
hold on;
    plot(0.001*heat_input, pipe_angle(distance_min_index2(19,:)),'k--', 'MarkerSize', 50);
hold on;
    plot(0.001*heat_input, pipe_angle(distance_min_index2(25,:)),'b-.', 'MarkerSize', 50);
hold on;
    plot(0.001*heat_input, pipe_angle(distance_min_index2(31,:)),'r', 'MarkerSize', 50);
hold on;

xlabel('Q (kW)');
xlim([0.5 50]);
ylim([3.5 4.5]);
ylabel('\phi*', 'fontweight','bold', 'fontsize', 16);
legend('d_0 = 70 mm','d_0 = 80 mm','d_0 = 90 mm','d_0 = 100 mm');

%% 

heat_input = linspace(500,50000,31);

[h_max, h_max_index] = max(h(:,:,:,:));
h_max = reshape(h_max, length(inlet_diameter), length(heat_input));
h_max_index = reshape(h_max_index, length(inlet_diameter), length(heat_input));

plot3(inlet_diameter, heat_input, h_max);
%% Other plots

[mind, mindi] = min(distance);
mind2d = reshape(mind, length(inlet_diameter), length(massflow)); % massflow?
mindi2d = reshape(mindi, length(inlet_diameter), length(massflow));

angles = pipe_angle(mindi2d);

figure;
surf(100.*inlet_diameter, massflow, angles');
% Why did I need to transpose this matrix to get the trends I "expected"?
% Are my optimizations completely off? Is heat input the one causing the 
xlabel('Inlet Diameter (cm)');
ylabel('Mass Flow Rate (kg/s)');
zlabel('Utopia Angle (degrees)');
set(gca,'FontSize',12)