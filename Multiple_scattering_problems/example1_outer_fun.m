% outer function to compute  amultiple scattering method which uses AG code
% for the left hand side matrix computation and PIM and midpoint for
% evaluating the right hand side.

clear all

% adding path to Andrews code
addpath('/Users/ojp18/Dropbox/PhD_DB/MATLAB/Andrew_Gibbs_new/HNABEMLAB-master/Examples/OP_functions')

% adding path to integrators etc
addpath('./General_functions')

% General set up

% Geometric set up
vertices1 = [-2*pi 0;
    0 0];
Gamma1=Screen(vertices1);

vertices2 = [200*pi 0;
    202*pi 0];
            
Gamma2=Screen(vertices2);

% General set up
%wavenumber
kwave=10;

theta = pi/4; %% need to change so it is consistent
d = [sin(theta) -cos(theta) ];
uinc=planeWave(kwave,d);

pMax = 4; %polynomial degree
cL = 2; %layers of grading per polynomial degree
sigmaGrad=0.15; %grading ratio
nLayers = cL*(pMax+1)-1; %number of layers of grading
throwAwayParam = 0; %no need to remove any basis elements
OverSample = 1.25; %choose amount to oversample by (40% here)

[v_N1, GOA1, colMatrix1, colRHS1, col_points1,...
v_N2, GOA2, colMatrix2, colRHS2, col_points2] ...
    = AG_code_pulling_out_info(pMax, cL, sigmaGrad, nLayers, OverSample, ...
    Gamma1, Gamma2, kwave, uinc );

% Info needed for our solve

% Isolating the collocation points in cartesian coordinates
L1 = sqrt( (vertices1(2, 1) - vertices1(1, 1))^2 + (vertices1(2, 2) - vertices1(1, 2))^2 );
x1_col = vertices1(1, 1) + col_points1*(vertices1(2, 1) - vertices1(1, 1))/L1;
y1_col = vertices1(1, 2) + col_points1*( vertices1(2, 2) - vertices1(1, 2) )/L1;

L2 = sqrt( (vertices2(2, 1) - vertices2(1, 1))^2 + (vertices2(2, 2) - vertices2(1, 2))^2 );
x2_col = vertices2(1, 1) + col_points2*(vertices2(2, 1) - vertices2(1, 1))/L1;
y2_col = vertices2(1, 2) + col_points2*( vertices2(2, 2) - vertices2(1, 2) )/L1;

G1 = [vertices1(1, 1), vertices1(1, 2), vertices1(2, 1), vertices1(2, 2)];
G2 = [vertices2(1, 1), vertices2(1, 2), vertices2(2, 1), vertices2(2, 2)];

% dui_dn1 = @(nodes) duidn(vertices1, L1, kwave, d, nodes);
% dui_dn2 = @(nodes) duidn(vertices2, L2, kwave, d, nodes);

N_approx = 2^(-6);
[x1, y1, t1, t1_mid, h1, h1vector, N1, L1] = discretisation_variables(G1, N_approx, k);
[x2, y2, t2, t2_mid, h2, h2vector, N2, L2] = discretisation_variables(G2, N_approx, k);


% Step 0


% Step 1

% Step 2
