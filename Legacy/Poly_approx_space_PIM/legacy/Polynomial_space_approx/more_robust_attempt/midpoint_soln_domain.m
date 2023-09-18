function z = midpoint_soln_domain(x1, x2, G, aj, k, Q)
% This is the midpoint function which will evaluate the integral \int_{Gj}
% \phik(x, y) \phi_{j}^{(?)}(y) ds (y). We are assuming that we have been
% given a point \xb = (x, y)^T and we are computing the value of the integral at
% that point
% 
% x is a point in the domain
% y is a point on the screen that we are integrating with respect to.
%
% Problem parameters:
% x1 is the x value of the point we are integrating
% x2 is the y value at the point we are integrating
% G is the coordinates of the screen
%
% Discretisation parameters:
% Q is the number of quadrature points used, quad discretisation

a = G(1);
b = G(2);
c = G(3);
d = G(4);
L = sqrt( (G(3) - G(1))^2 + (G(4) - G(2))^2 );  % length of G1
% h = L/N;

hq = Q/L;  % quadrature step size and weight as midpoint rule

t = [hq/2:hq:L - hq/2]; % integration nodes
z = 0;
for j = 1:length(t)
    
    z = z + hq*1i*besselh(0,k *sqrt( ( x1 - a - t*(c-a)/L )^2 ...
    + ( x2 - b - t*(d - b)/L )^2 ) )*coeff_2_soln_midpoint_individual(aj, L, t(j), N)/4;
    
end
