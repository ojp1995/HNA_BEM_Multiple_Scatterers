function [us, x1, y1] = compute_scattered_field_beam(k, X, Y, nqx, nqy,...
    h, fnq, G, nq_1D)
% In this function we will compute the scattered field from a beam source.
% This could be an incident source or it could be the scattered field that
% we then minus the plane wave from to find the total field
%
% Needs work in this section.

us = zeros(length(Y), length(X));

for ix = 1:length(X)
    
    for iy = 1:length(Y)
        
        dist = sqrt( (X(ix) - nqx).^2 + (Y(iy) - nqy).^2);
        
        us(iy, ix) = 1i*h*sum(besselh(0, k*dist).*fnq.')/4;
        
    end
    
end

L = sqrt( (G(3) - G(1))^2 + (G(4) - G(2))^2 );
x1 = G(1) + nq_1D*(G(3) - G(1))/L;
y1 = G(2) + nq_1D*( G(4) - G(2) )/L;