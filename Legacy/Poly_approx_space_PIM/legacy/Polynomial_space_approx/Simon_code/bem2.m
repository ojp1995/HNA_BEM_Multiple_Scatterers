function [dudn, b] = bem2(x,y,h,k,d)
%
% function dudn = bem2(x,y,h,k,d)
%
% N.B. function bem carries out the same calculations, but slightly less
% efficiently (it is less vectorised).
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
% This function, together with the function us, solve this scattering
% problem by a simple boundary element method. It is assumed that Gamma is
% divided (exactly or approximately) into a number of straight-line
% boundary element methods. The function us computes the scattered field by
% a discretisation of the representation formula
%
% u(rr) = u^i(rr) - \int_Gamma \Phi(rr,rr_s) dudn(rr_s) ds(rr_s),    (1)
%
% where dudn is the normal derivative (the normal directed into D) of u on
% Gamma, and
%
% Phi(rr,rr_s) := (i/4}H_0^{(1)}(k|rr-rr_s|)
% 
% is the standard fundamental solution of the Helmholtz equation. This
% function bem2 computes values of dudn numerically by a simple piecewise
% constant boundary element method, collocating (1) at the midpoint on Gamma 
% of each element (where u = 0).
%
% The inputs are:
%
% x,y        real vectors of length some lenth N
%            (x(j),y(j)) are the coordinates of the midpoint of the jth
%            element.
%
% h          real vector of length N
%            h(j) is the length of the jth element.
%
% k          real or complex scalar, with Re k > 0, Ik k >= 0.
%            k is the wave number.
%
% d          real vector of length 2.
%            d is a unit vector in the direction of the plane incident
%            wave.
%
% The output is:
%
% dudn       complex vector of length N.
%            dudn(j) is an approximation to the value of the normal derivative, dudn
%            at the midpoint of the jth element.
%
N = length(x); % Number of boundary elements.
A = zeros(N,N); b = zeros(N,1); % A and b will be the matrix and RHS in the discretised boundary integral equation.
Nsing = 200; % Number of quadrature points to use when calculating the diagonal elements of A.
gr = (1:Nsing)-0.5; 

b = ui(x,y,k,d(1),d(2)).'; % Evaluating the RHS of the BIE (1) at the collocation points om Gamma.
for n = 1:N
    % Evaluating the nth row of A
    h_small = 0.5*h(n)/Nsing; % Stepsize for composite midpoint rule evaluation of the diagonal element
    A(n,n) = (i/2)*h_small*sum(besselh(0,k*gr*h_small)); % First the entry on the diagonal  
    select = (1:N)~=n; % Now the off-diagonal elements
    A(select,n) = (i*h(n)/4)*besselh(0, k*sqrt( (x(select)-x(n)).^2 + (y(select)-y(n)).^2 ) );
end
dudn = A\b; 