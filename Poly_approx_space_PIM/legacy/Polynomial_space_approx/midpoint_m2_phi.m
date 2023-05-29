function I = midpoint_m2_phi(a, b, x, FUN_m2, FUN_phi, N)
% This is the midpoint function which will evaluate the integral of FUN(x)
% between a and b (where support is)
%
% Problem parameters:
% a is start of the integral
% b is the end of the integral
% FUN_m2 is the function m2(s) = k(s) - sigma1(s)m1_tilde(s)
% FUN_phi is the basis function, piecewise constants
% x is the fixed point we are evlauating the integral at
%
% Discretisation parameters:
% N is the numbe rof quadrature points we are evaluating this integral with
%

Q = N+1; %% got around bodge by removing lines below, think this is ok, if not change back
%%%% BODGE!!! Have had to increase the number of intervals, ie Q = N+1 to
%%%% prevent a singularity

h = (b - a)/Q;  % discretisation of interval [a, b]
n =  [a + h/2 : h: b - h/2];  % midpoint quadrature nodes

I = 0;  % Initialisation of summing variable for integration


for j = 1:Q
    I = I + h*FUN_m2(abs(x - n(j))); %%*FUN_phi(a, b, n(j)); - REMOVED as not needed and doesn't make sense
    
%     if FUN_phi(a, b, n(j)) ==0
%         keyboard
%     else
%     end
end
% z = sum(w*FUN_m2(x - n)*FUN_phi(a, b, n));