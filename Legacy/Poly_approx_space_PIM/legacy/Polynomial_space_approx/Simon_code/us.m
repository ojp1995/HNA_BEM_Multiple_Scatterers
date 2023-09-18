function u = us(X,Y,k,x,y,h,dudn)
%
% function u = us(X,Y,k,x,y,h,dudn)
%
% Let D be the connected (but not necessarily simply-connected) exterior of
% a bounded Lipschitz domain in 2D, and let Gamma be its boundary. Consider
% the following time harmonic acoustic scattering problem: given an
% incident plane wave u^i(rr) = exp(i k rr.d), where k is the wave number, rr a
% position vector, and d is a unit vector in the direction of the plane wave, find the total
% wave u that satisfies the Helmholtz equation
%
% Delta u + k^2 u = 0           in D,
%
% u = 0                         on Gamma,
%
% and the Sommerfeld radiation condition
%
% \partial u^s/\partial r = i k u^s + o(r^(-1/2))
%
% as r = |rr| \to \infty, where u^s := u - u^i is the scattered field.
%
% This function, together with the function bem2, solve this scattering
% problem by a simple boundary element method. It is assumed that Gamma is
% divided (exactly or approximately) into a number of straight-line
% boundary element methods. This function us computes the scattered field by
% a discretisation of the representation formula
%
% u(rr) = u^i(rr) - \int_Gamma \Phi(rr,rr_s) dudn(rr_s) ds(rr_s),    (1)
%
% where dudn is the normal derivative (the normal directed into D) of u on
% Gamma, and
%
% Phi(rr,rr_s) := (i/4}H_0^{(1)}(k|rr-rr_s|)
% 
% is the standard fundamental solution of the Helmholtz equation. The
% function bem2 computes values of dudn numerically by a simple piecewise
% constant boundary element method, collocating (1) at the midpoint on Gamma 
% of each element (where u = 0).
%
% The inputs are:
%
% X, Y       real scalars, vectors, or matrices (with X and Y having the same dimensions)
%            These inputs contain the (X,Y) coordinates of a single point,
%            or of a vector or matrix of points at which one wishes to calculate
%            u^s.
%
% k          real or complex scalar, with Re k > 0, Ik k >= 0.
%            k is the wave number.
%
% x,y        real vectors of length some lenth N
%            (x(j),y(j)) are the coordinates of the midpoint of the jth
%            element.
%
% h          real vector of length N
%            h(j) is the length of the jth element.
%
% dudn       complex vector of length N.
%            dudn(j) is an approximation to the value of the normal derivative, dudn
%            at the midpoint of the jth element.
%
% The output is:
%
% u          complex array of same size as X and Y
%            Contains numerical approximations to the values of the
%            scatteted field u^s at the (X,Y) coordinates stored in XX and
%            YY.
%

u = zeros(size(X));
for m = 1:length(h)
    dist = sqrt( (X-x(m)).^2 + (Y-y(m)).^2 );
    u = u - (i*h(m)*dudn(m)/4)*besselh(0,k*dist);
end