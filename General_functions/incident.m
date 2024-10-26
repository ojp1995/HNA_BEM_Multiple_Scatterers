function z = incident(k, theta, x1, y1)
% robust adapting method
%
% here we will compute the plane wave incident for a screen with
% coordinates G.

% Problem parameters:
% x = (x1, y1), the corrodinates of the screen
% k is the wavenumber
% theta is the angle between the downwards vertical and the incident wave anticlockwise
% G is the coordinates of the screen
% x are the points we are evaluating this at
%
%
% Discretisation parameters:
% N is the discretisations

% L = sqrt( (G(3) - G(1))^2 + (G(4) - G(2))^2 );  % length of the screem
% h = L/N;



% x1 = G(1) + x.*(G(3) - G(1))/L;
% x2 = G(2) + x.*( G(4) - G(2) )/L;

z = exp( 1i*k*( x1*sin(theta) - y1*cos(theta) ) );