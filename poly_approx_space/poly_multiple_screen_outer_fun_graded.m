% outer function for the polynomial solver with graded quadrature
clear all

% only adding one path to centralise solvers
addpath('../General_functions/')

% switch and things get interesting
G1 = [-2*pi, 2*pi, 0, 0];
L1 = sqrt( (G1(3) - G1(1))^2 +(G1(4) - G1(2))^2 );

G2 = [2*pi, 0, 5*pi, 3*pi]; 
L2 = sqrt( (G2(3) - G2(1))^2 +(G2(4) - G2(2))^2 );

C_wl= 2^-6;

k = 10;  % wavenumber

theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

% now compute the 4 matrices (hopefully with a basis function support 
% being general)

% grid for support of basis functions
% This will almost all be changed so are very much holding lines to make
% sure that it works for the most simple case before we try anything more
% complicated
h1 = L1/1000;
h2 = L2/1000;
t1_grid = [0: h1: L1];
t2_grid = [0: h2: L2];

[S11, S12, S21, S22, u_inc1, uinc2] = ...
    compute_matrices_for_iterative_solve(G1, G2, k, )


