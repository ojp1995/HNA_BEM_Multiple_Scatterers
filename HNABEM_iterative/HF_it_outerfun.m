% Iterative method in a HNA space
% Paths to general functions and BEAM_HNABEMLAB
addpath('../General_functions/')
addpath('../../BEAM_HNABEMLAB/')
addPathsHNA  % allows HNABEM to find all of the relevatn subfolders

%% Geometrical set up specific to HNABEMLAB set up

vertices1 = [-2*pi 2*pi;
    0, 0];

vertices2 = [2*pi 0;
    5*pi 3*pi];

Gamma1=Screen(vertices1);
Gamma2=Screen(vertices2);

kwave=10;
theta = 0;
d = [sin(theta) -cos(theta) ];
uinc=planeWave(kwave,d);

pMax = 4; %polynomial degree
cL = 2; %layers of grading per polynomial degree
sigmaGrad=0.15; %grading ratio
nLayers = cL*(pMax+1)-1; %number of layers of grading
throwAwayParam = 0; %no need to remove any basis elements
OverSample = 1.25; %choose amount to oversample by (40% here)

[v_N1, GOA1, colMatrix1, colRHS1, col_points1,...
v_N2, GOA2, colMatrix2, colRHS2, col_points2, VHNA1, VHNA2] ...
    = AG_code_pulling_out_info(pMax, cL, sigmaGrad, nLayers, OverSample, ...
    Gamma1, Gamma2, kwave, uinc );

%% Info needed for our part of the solve

