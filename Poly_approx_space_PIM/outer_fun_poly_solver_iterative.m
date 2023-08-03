% outer function for polynomial solver using as many functions from HNA
% folder as possiible

clear all

% adding path to solvers in HNA folder
addpath()

% switch and things get interesting
G1 = [-2*pi, 2*pi, 0, 0];

G2 = [2*pi, 0,  5*pi, 3*pi]; 

C_wl= 2^-6;

k = 10;  % wavenumber

theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

% number of iterations for approximation and truth

R_it = 20;  % R_it = 13 (I think) for len =0.5. R_it > 50, len = 2;
R_true = 20;

