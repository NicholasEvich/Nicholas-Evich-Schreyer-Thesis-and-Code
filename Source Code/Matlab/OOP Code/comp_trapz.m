% Nicholas Evich
% 2/11/2020
% MATH 456
% Numeric integration function utilizing the composite trapezoid rule

function [int] = comp_trapz(f, a, b, N)
% Composite trapezoid rule function for numeric intergration
    % f = function handle f(x)
    % a = starting point
    % b = end point
    % N = number of panels
h = (b - a)/N;

int = 0.5*(f(a) + f(b));

for i = 1:N-1
   int = int + f(a + i*h);
end

int = int*h;

end

