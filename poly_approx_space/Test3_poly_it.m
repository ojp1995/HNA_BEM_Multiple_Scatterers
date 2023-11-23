% test 3 - partial shadowing

clear all

addpath('../General_functions/')  % access to solvers needed

% introducing the screens, storing the data in a struct object 

G1_data.G = [0, 0, 6*pi, 0];

G2_data.G = [1, 1, 1+ 6*pi, 1]; 

% coefficients needed for creating grid for basis functions and quadrature
% points
Lgrad_coeff = 0.15;
alpha = 2;

% creating basis function information
C_wl_bf1 = 1/10;
C_wl_bf2 = 1/10;

C_wl_quad= 1/20;

R_max = 20;

k = 10;

theta = ;
C1 = 1;
C2 = pi;


% solve
[G1_data, G2_data, aj_1_R, aj_2_R, us] = ...
    compute_iteratuve_poly_scattering_prob_2_screens(G1_data, G2_data, ...
    k, Lgrad_coeff, alpha, C_wl_bf1, C_wl_bf2, C_wl_quad, R_max, theta, ...
    C1, C2, true, true);
