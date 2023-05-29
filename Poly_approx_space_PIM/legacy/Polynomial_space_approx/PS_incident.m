function z = PS_incident(G, L, k, s, x, y )

% Problem paramters:
% G is the coordinates of the screen:
% (a, b) is the start point of the screen
% (c, d) is the end point of the screen
% 
% L is the length of the screen
% s is the vector of collocation points
% k is the wavenumber 
% (x, y) is the point source
a = G(1);
b = G(2);
c = G(3);
d = G(4);

x1 = a + (s(1)./L).*(c-a) ; y1 = b + (s(1)./L).*(d-b);

%%% WHAT IS THIS A FOR???
A=exp(-1i*k*y1);

z = A*besselh(0, k*sqrt( (a + (s./L).*(c-a)  - x).^2 + ( b + (s./L).*(d-b) - y ).^2 ) );