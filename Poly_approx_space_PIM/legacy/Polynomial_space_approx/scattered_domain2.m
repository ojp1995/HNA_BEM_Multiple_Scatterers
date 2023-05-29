function [us, x1, y1, x2, y2] = scattered_domain2( G1, G2, aj1, aj2, k, X, Y, N1, N2)
% In this function we will compute the scattered domain for any given
% coefficients on the same given screen

% problem parameters
% G are the coordinates of the screen
% aj are the coefficients for the solution we are interested in
% k is the wavenumber
% discretisation parameters
% X is the discretisation in the X direction
% Y is the discretisation in the Y direction
% N is the number of intervals on the screen

us = zeros(length(Y), length(X));
L1 = sqrt( (G1(3) - G1(1))^2 + (G1(4) - G1(2))^2 );
h1 = L1/N1;  % computing the midpoint weights, discretisation interval.
t1 = [h1/2:h1:L1-h1/2];  % midpoints, quad nodes

L2 = sqrt( (G2(3) - G2(1))^2 + (G2(4) - G2(2))^2 );
h2 = L2/N2;  % computing the midpoint weights, discretisation interval.
t2 = [h2/2:h2:L2-h2/2];  % midpoints, quad nodes

for ix = 1:length(X)% loop through every point in x direction
    
    for iy = 1:length(Y)  % looops through every point in the y direction
        
        us(iy, ix) = aj1.'*midpoint_domain_int(k, G1(1), G1(2), G1(3), G1(4), X(ix), Y(iy), L1, h1, t1).' ...
            + aj2.'*midpoint_domain_int(k, G2(1), G2(2), G2(3), G2(4), X(ix), Y(iy), L2, h2, t2).';
        
%         us(iy, ix) = midpoint_soln_domain(X(ix), Y(iy), G, aj, k, Q);
        
    end
    
end

x1 = G1(1) + t1*(G1(3) - G1(1))/L1;
y1 = G1(2) + t1*( G1(4) - G1(2) )/L1;

x2 = G2(1) + t2*(G2(3) - G2(1))/L2;
y2 = G2(2) + t2*( G2(4) - G2(2) )/L2;

