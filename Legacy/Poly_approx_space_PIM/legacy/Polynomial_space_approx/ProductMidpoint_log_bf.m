function I = ProductMidpoint_log_bf(a, b, k, x, m, phi, N)
% We are integrating I(x) = \int_{a}^{b} m(x - t) ln|x - t| \phi(t) dt using the
% product integration rule with the midpoint method, we are approximating
% the integral I(x) as I(x) \approx \sum_{j = 1}^{N} m(x - n_{j}) \phi(n_j)
% \int_{x_{j}}^{x_{j+1}} ln(x - t) d t, for some x either in the interval
% or not.
% 
% Have added in a line for the support of the basis functions
%
% Problem parameters:
% a is the lower bound of the integral
% b is the upper bound of the integral
% x is a fixed point
% m is the smooth function, m(x, t)
% phi is the function phi(t)
%
% Discretisation paraemters:
% N is the number of intervals the integral is being split up into
%
% Assumption:
% 1. x is a scalar
% 2. may need to look at a different formula for when x is far away from
% [a, b].


h = (b - a)/N;  % discretisation of interval [a, b]
n =  [a + h/2 : h: b - h/2];  % midpoint quadrature nodes

I = 0;  % Initialisation of summing variable for integration

for j = 1:N
    
    I = I + m(abs(x - n(j)))*phi(a + (j-1)*h, a+j*h, n(j))*v1_weights_mid(a + (j-1)*h, a+j*h, k, x);
    
end









