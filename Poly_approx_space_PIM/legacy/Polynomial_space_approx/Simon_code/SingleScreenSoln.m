%
% Computation by a piecewise constant collocation method of the solution 
% phi on a straight line screen Gamma when
%
% S phi = u^i|_Gamma,
% 
% and S is the single-layer potential integral operator, with kernel
%
% Phi(r,r_0) = (i/4) H_0^{(1)}(k|r-r_0|)
%
% and u_i is the incident plane wave
%
format compact
k = 1 % the wavenumber
a = 0; b = 2*pi; % the screen is between (a,0) and (b,0)
d = [0.7071,-0.7071]; % the direction of the incident wave
N = 100 % the number of boundary elements
h = (b-a)/N % the length of each element
x = h*((1:N)-0.5); y = zeros(size(x)); % the x and y coordinates of the element midpoints
hh = h*ones(size(x)); % the lengths of the elements
[phi, A] = bem2(x,y,hh,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
 
figure(); plot(x, real(phi));
title('Simon code, single screen')
xlabel('x')
ylabel('\phi(x)')