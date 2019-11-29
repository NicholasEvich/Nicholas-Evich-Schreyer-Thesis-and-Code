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