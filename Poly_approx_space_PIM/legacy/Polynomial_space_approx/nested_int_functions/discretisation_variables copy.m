function [x, y, t, h, hvector, N, L] = discretisation_variables(G, C_wl, k)
% In this function we will produce all the discretisation paramters needed
% for the problem, given the following information:
% G1, the coordinates for the screen
% C_wl, the number of degrees of freedom per wavelength
% k the wavenumber

% first computing the length of the screen
L = sqrt( (G(3) - G(1))^2 + (G(4) - G(2))^2 ); 

% number of dof
N = ceil(k*L./(C_wl*2*pi));

% step size
h = L/N;

% computing the midpoints coordinates for \Gamma
x = G(1)+((1:N)-0.5)*(G(3)-G(1))/N;
y = G(2)+((1:N)-0.5)*(G(4)-G(2))/N;

%paramterisation of the screen
t = [0:h:L];

% step size of screen
hvector = h*ones(size(x)); % Lengths of the elements on screen 1 









