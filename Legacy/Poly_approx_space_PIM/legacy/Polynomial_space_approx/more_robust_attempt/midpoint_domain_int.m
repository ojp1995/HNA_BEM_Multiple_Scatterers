function z = midpoint_domain_int(k, a, b, c, d, x1, x2, L, w, t)
% new function to try to integrate
% \int_{t_{j-1}}^{t_{j}} H_{0}^{(1)}(k \vert \xb - \yb \vert) ds (\yb)
%
% Problem parameters:
%  k is the wavenumber
% (a, b) is the start coordinates of the screen
% (c, d) is the end of the screen
% (x1, x2) are the coordinates the points in the domain
% L is the length of the screen
% 
% Discretisation parameters:
% w are the weights, in this case h
% t are the intergation nodes

z = 1i*w*besselh(0, k*( sqrt( ( x1 - a - t*(c - a)/L ).^2 ...
    + ( x2 - b - t*(d - b)/L ).^2) ) )/4;