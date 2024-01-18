function z = duidn(G_data, L, k, d, t)
% In this function we will be computing the the normal derivative of the
% incident wave of the screen due to a plane wave incident on it, where the
% plane wave has vector direction d. We will be evaluating this at the
% point x = (x1, y1).
%
% Problem parameters:
% Gamma - The coordinates of the screen, (a, b ; c, d)
% L is the length of the screen
% k the wavenumber
% d is the direction of the incident wave.
% t is the vector of points we are evaluating this at, we are assuming that
% it has been parameterised before hand

% Caution, not sure if this can handle vectors for the x yet.

% x1 = G(1, 1) + t.*(G(2, 1) - G(1, 1))/L;
% x2 = G(1, 2) + t.*( G(2, 2) - G(1, 2) )/L;

x1 = G_data.G(1) + t.*(G_data.G(3) - G_data.G(1))/G_data.L;
x2 = G_data.G(2) + t.*(G_data.G(4) - G_data.G(2))/G_data.L;



% think there is an error below, should be + on the outside of the term
% then - infront of d(1)*G(2, 2) ... 
% z = 1i*k*exp(1i*k*( x1*d(1) + x2*d(2)) )*( d(1)*(G(2, 2) - G(1, 2)) ...
%     + d(2)*(G(2, 1) - G(1, 1)) )/L;

z = 1i*k*( (d(1)*G_data.n(1)*exp(1i*k*x1*d(1))) ...
    + (d(2)*G_data.n(2)*exp(1i*k*x2*d(2))) );

