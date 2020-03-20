% Nicholas Evich
% 2/11/2020
% MATH 456
% Function for computing a degree-N least squares approximation of the
% function f with user defined weighting functions, numeric integration
% function, and interval over which to evaluate

function [A] = least_square_poly(f, w, N, a, b, int, panels)
% Least-squares polynomial approximation function
    % f = function to be approximated f(x)
    % w = weighting function
    % N = desired degree of least square polynomial
    % a = start point
    % b = end point
    % int = numeric integration function
    % panels = number of steps to use in numeric integration strategy
    
    A = zeros(1, N+1); % empty vector of approximation coefficients
    
    phi_k2 = @(x) 1;
    integ_num = @(x) (x*w(x)*(phi_k2(x))^2);
    integ_den = @(x) (w(x)*(phi_k2(x))^2);
    
    B1 = int(integ_num, a, b, panels)/int(integ_den, a, b, panels);
    
    phi_k1 = @(x) x - B1;
    
    % Calculating A0
    integ_num = @(x) w(x)*phi_k2(x)*f(x);
    integ_den = @(x) w(x)*(phi_k2(x))^2;
    A(1) = int(integ_num, a, b, panels)/int(integ_den, a, b, panels);
    
    % Calculating A1
    integ_num = @(x) w(x)*phi_k1(x)*f(x);
    integ_den = @(x) w(x)*(phi_k1(x))^2;
    A(2) = int(integ_num, a, b, panels)/int(integ_den, a, b, panels);
    
    if N >=2
       for k = 3:N+1 % start at 3 because matlab array indexing is stupid
           integ_num = @(x) x*w(x)*(phi_k1(x))^2;
           integ_den = @(x) w(x)*(phi_k1(x))^2;
           Bk = int(integ_num, a, b, panels)/int(integ_den, a, b, panels);
           integ_num = @(x) x*w(x)*phi_k1(x)*phi_k2(x);
           integ_den = @(x) w(x)*(phi_k2(x))^2;
           Ck = int(integ_num, a, b, panels)/int(integ_den, a, b, panels);
           phi_k = @(x) (x - Bk)*phi_k1(x) - Ck*phi_k2(x);
           
           integ_num = @(x) w(x)*phi_k(x)*f(x);
           integ_den = @(x) w(x)*(phi_k(x))^2;
           A(k) = int(integ_num, a, b, panels)/int(integ_den, a, b, panels);
           
           % Setting correct phi functions for the next iteration
           phi_k2 = @(x) phi_k1(x);
           phi_k1 = @(x) phi_k(x);
       end
    end
    
end

