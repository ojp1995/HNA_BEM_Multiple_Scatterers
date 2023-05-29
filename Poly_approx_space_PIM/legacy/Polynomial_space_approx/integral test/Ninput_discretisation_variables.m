function [x, y, t, t_mid, h, hvector, L] = Ninput_discretisation_variables(G, N, k)
% In this function we will produce all the discretisation paramters needed
% for the problem, given the following information:
% G1, the coordinates for the screen
% N, the number of degrees of freedom 
% k the wavenumber

% first computing the length of the screen
L = sqrt( (G(3) - G(1))^2 + (G(4) - G(2))^2 ); 

% step size
h = L/N;

% computing the midpoints coordinates for \Gamma
x = G(1)+((1:N)-0.5)*(G(3)-G(1))/N;
y = G(2)+((1:N)-0.5)*(G(4)-G(2))/N;

%paramterisation of the screen
t = [0:h:L];

% midpoints
t_mid = [ h/2: h: L - h/2];

% step size of screen
hvector = h*ones(size(x)); % Lengths of the elements on screen 1 









