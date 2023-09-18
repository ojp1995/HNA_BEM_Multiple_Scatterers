function [us, x1, y1] = scattered_domain( G, aj, k, X, Y, N)
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
L = sqrt( (G(3) - G(1))^2 + (G(4) - G(2))^2 );
h = L/N;  % computing the midpoint weights, discretisation interval.
t = [h/2:h:L-h/2];  % midpoints, quad nodes

for ix = 1:length(X)% loop through every point in x direction
    
    for iy = 1:length(Y)  % looops through every point in the y direction
        
        us(iy, ix) = aj.'*midpoint_domain_int(k, G(1), G(2), G(3), G(4), X(ix), Y(iy), L, h, t).';
        
%         us(iy, ix) = midpoint_soln_domain(X(ix), Y(iy), G, aj, k, Q);
        
    end
    
end

x1 = G(1) + t*(G(3) - G(1))/L;
y1 = G(2) + t*( G(4) - G(2) )/L;