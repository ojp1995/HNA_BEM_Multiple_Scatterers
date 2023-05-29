function z = PS_incident(k, xs, ys, tx, ty)
% this function will compute the point source starting from an origin (xs,
% ys) incident on a screen \Gamma with wavenumber k.
%
% Problem parameters
% (tx, ty) vector of parameterised points along the screen
% (xs, ys) is the origin of the point source
% k is the wavenumber.

 z = besselh(0, k*sqrt( (tx - xs).^2  + (ty - ys).^2));
