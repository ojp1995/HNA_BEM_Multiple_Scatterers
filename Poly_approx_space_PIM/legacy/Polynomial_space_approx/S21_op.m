function [S21, col_points] = S21_op(G1, G2, k, N1, N2)
% This function computes S21\phi_{1}^{(r)} = \int_{\G1}\phik(x,
% y)\phi_{1}^{(r)}(y) ds(y), x \in \G2
% Idea here is the beam originating from G1 incident on G2.
%
% Problem parameters
% G1 is the coordinates of G1
% G2 are the coordinates of G2
% k is the wavenumber
%
% Discretisation parameters:
% t is the integration parameter
% s is the collocation point
% h is the intergration weights
% N1 is the number of intervals on G1
% N2 is the number of intervals on G2

% coordinates for screens
a1 = G1(1);
b1 = G1(2);
c1 = G1(3);
d1 = G1(4);
L1 = sqrt( (G1(3) - G1(1))^2 + (G1(4) - G1(2))^2 );  % length of G1
h1 = L1/N1;

a2 = G2(1);
b2 = G2(2);
c2 = G2(3);
d2 = G2(4);
L2 = sqrt( (G2(3) - G2(1))^2 + (G2(4) - G2(2))^2 );  % length of G2
h2 = L2/N2;

% computing the t parameter, these will be the midpoints of the intervals
% on G1
t = [h1/2: h1: L1 - h1/2];

% computing the s parameter, the midpoints of the intervals on G2
col_points = [h2/2: h2: L2 - h2/2];

for j = 1:N1  % basis function loop, t param
    
    for l = 1:length(col_points)  % collocation loop, s param
        
        S21(l, j) = h1*1i*besselh(0, k*sqrt( (a2 + col_points(l)*(c2 - a2)/L2  - ( a1 + t(j)*(c1 - a1)/L1 ))^2 ...
            +  (b2 + col_points(l)*(d2 - b2)/L2 - (b1 + t(j)*(d1 - b1)/L1) )^2 ))/4;
        
    end
    
end