function [x, y, t, t_mid, h, hvector, N, L] = discretisation_variables_edge(G, C_wl, k)
%
%
% Adapted so it is away from the edges avoid singularity problems
%
% In this function we will produce all the discretisation paramters needed
% for the problem, given the following information:
% G1, the coordinates for the screen
% C_wl, the number of degrees of freedom per wavelength
% k the wavenumber
tol = k/10000;
% first computing the length of the screen
G(1) = G(1) + tol;
G(2) = G(2) - tol;
G(3) = G(3) - tol;
G(4) = G(4) + tol;

L = sqrt( (G(3) - G(1))^2 + (G(4) - G(2))^2 ); 

% number of dof
N = ceil(k*L./(C_wl*2*pi));

% step size
h = L/N;

tol = k*h/2;  % tolerance we are moving away from the end points

% computing the midpoints coordinates for \Gamma
x = (G(1))+((1:N)-0.5)*((G(3)) - (G(1)))/N;
y = (G(2))+((1:N)-0.5)*((G(4)) - (G(2)))/N;

%paramterisation of the screen
t = [tol:h:L - tol];

t_mid = [tol + h/2: h: L - h/2 - tol];

% step size of screen
hvector = h*ones(size(x)); % Lengths of the elements on screen 1 









