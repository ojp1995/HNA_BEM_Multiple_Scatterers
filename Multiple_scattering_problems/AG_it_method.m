% In this function we will be using Andrews code entirely to perform the
% iterative method.

clc;
clear classes;
close all;
clear all;
% addpath ../HNABEMLAB;
addpath('/Users/ojp18/Dropbox/Mac/Documents/GitHub/BEAM_HNABEMLAB')
addpath('/Users/ojp18/Dropbox/Mac/Documents/GitHub/HNA_BEM_Multiple_Scatterers/General_functions')
addpath('/Users/ojp18/Dropbox/Mac/Documents/GitHub/BEAM_HNABEMLAB/Beam/Beam source stuff')
addPathsHNA;

%================== whatev
kwave = 10.0;

%create 'screen' object ---------------------------------------------------
vertices1 = [-2*pi 2*pi;
    0, 0];
Gamma1=Screen(vertices1);

vertices2 = [2*pi 0;
    5*pi 3*pi];
            
Gamma2=Screen(vertices2);

% General set up
%wavenumber
kwave=10;
theta = 0;
% theta = pi/4; %% need to change so it is consistent
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

% Isolating the collocation points in cartesian coordinates
L1 = sqrt( (vertices1(2, 1) - vertices1(1, 1))^2 + (vertices1(2, 2) - vertices1(1, 2))^2 );
x1_col = vertices1(1, 1) + col_points1*(vertices1(2, 1) - vertices1(1, 1))/L1;
y1_col = vertices1(1, 2) + col_points1*( vertices1(2, 2) - vertices1(1, 2) )/L1;

L2 = sqrt( (vertices2(2, 1) - vertices2(1, 1))^2 + (vertices2(2, 2) - vertices2(1, 2))^2 );
x2_col = vertices2(1, 1) + col_points2*(vertices2(2, 1) - vertices2(1, 1))/L1;
y2_col = vertices2(1, 2) + col_points2*( vertices2(2, 2) - vertices2(1, 2) )/L1;

G1 = [vertices1(1, 1), vertices1(1, 2), vertices1(2, 1), vertices1(2, 2)];
G2 = [vertices2(1, 1), vertices2(1, 2), vertices2(2, 1), vertices2(2, 2)];

n1 = [-(G1(4) - G1(2)), G1(3) - G1(1)]/L1;
n2 = [-(G2(4) - G2(2)), G2(3) - G2(1)]/L2;

% dui_dn1 = @(nodes) duidn(vertices1, L1, kwave, d, nodes);
% dui_dn2 = @(nodes) duidn(vertices2, L2, kwave, d, nodes);

N_approx = 2^(-10);
[x1, y1, t1, t1_mid, h1, h1vector, N1, L1] = discretisation_variables(G1, N_approx, kwave);
[x2, y2, t2, t2_mid, h2, h2vector, N2, L2] = discretisation_variables(G2, N_approx, kwave);

% Infor needed for plotting
x1_plot_1D = [t1_mid(1): 0.001: t1_mid(end)];
x1_plot = linspace(x1_col(1), x1_col(end), length(x1_plot_1D)).';
y1_plot = linspace(y1_col(1), y1_col(end), length(x1_plot_1D)).'; 

away_edge = 3;
x2_plot_1D = [t2_mid(away_edge): 0.001: t2_mid(end - away_edge)];

% x2_plot_1D = linspace(t2_mid(1), t2_mid(1), 1e4).';
x2_plot = linspace(x2_col(away_edge), x2_col(end - away_edge), length(x2_plot_1D)).';
y2_plot = linspace(y2_col(away_edge), y2_col(end - away_edge), length(x2_plot_1D)).';



weights = v_N1.edgeComponent.nodes;
param_vals = v_N1.edgeComponent.weights;
density_samples = v_N1.eval(param_vals, 1) + GOA1.eval(param_vals, 1);

beam_inc_S21_phi1_0 = beamSol(kwave, density_samples, weights, param_vals, Gamma1);

S2=singleLayer(kwave,Gamma2);
[v_N_beam1, GOA_beam1, colMatrix_beam1, colRHS_beam1, T_beam1] = ColHNA(S2, VHNA2, beam_inc_S21_phi1_0.', Gamma2, 'oversample', OverSample, 'progress');

f_2_1_RHS = incident(kwave, theta, x2_col, y2_col) ...
    - beam_inc_S21_phi1_0.eval(x2_col, y2_col) - GOA2.eval(col_points2, 1)...
    - GOA_beam1.eval(col_points2, 1);


v_N_G2_r1 = compute_coeffs_given_A_and_f(colMatrix2,f_2_1_RHS , VHNA2);

phi2_1_manual = v_N_G2_r1.eval(x2_plot_1D.', 1) - GOA2.eval(x2_plot_1D.', 1)...
    - GOA_beam1.eval(x2_plot_1D.', 1);
%% Plotting different parts of the solution.
figure()
plot(x2_plot_1D, real(v_N2.eval(x2_plot_1D.', 1) ...
    + GOA2.eval(x2_plot_1D.', 1)))
title('Plot of solution on $\Gamma_{2}$ just due to plane wave incident')

figure()
plot(x2_plot_1D, real(v_N_beam1.eval(x2_plot_1D.', 1) ...
    + GOA_beam1.eval(x2_plot_1D.', 1)));
title('Plot of beam source, $S_{21}\phi_{1}^{0}$')

figure()
plot(x2_plot_1D, real(v_N_beam1.eval(x2_plot_1D.', 1)));
title('Plot of beam source, $S_{21}\phi_{1}^{0}$, Oscil part only')

figure()
plot(x2_plot_1D, real(GOA_beam1.eval(x2_plot_1D.', 1)));
title('Plot of beam source, $S_{21}\phi_{1}^{0}$, Leading order part only')


% combine the solutions
phi2_1 = v_N2.eval(x2_plot_1D.', 1) + GOA2.eval(x2_plot_1D.', 1) ...
    - (v_N_beam1.eval(x2_plot_1D.', 1) + GOA_beam1.eval(x2_plot_1D.', 1));

figure()
plot(x2_plot_1D, real(phi2_1));
title('Solution from single screen of uinc - beam source')

figure(); domainPlot(Gamma2, beam_inc_S21_phi1_0, GOA_beam1, v_N_beam1, kwave, 10000);
title('Plot in domain of $S_{21}\phi_{1}^{(0)}$ incident on $\Gamma_{2}$')
figure(); domainPlot(Gamma2, uinc, GOA2, v_N2, kwave, 10000);
title('Plot in domain of plane wave incident on $\Gamma_{2}$')


