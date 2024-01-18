function [G1_data, G2_data, phi1_r, phi2_r, v_N1cell, v_N2cell] = ...
    HF_it_outer_function(kwave, vertices1, vertices2, R_max, theta, ...
    C_wl_quad_outer, C_wl_quad_inner)
% The outer function to computer the iterative method given minimal inputs

% converting so suitable for Andrews solver
Gamma1=Screen(vertices1);
Gamma2=Screen(vertices2);

% infor for incident wave
d = [sin(theta) -cos(theta) ];
uinc=planeWave(kwave,d);

% variables needed for Andrews solver
pMax = 4; %polynomial degree
cL = 2; %layers of grading per polynomial degree
sigmaGrad=0.15; %grading ratio
nLayers = cL*(pMax+1)-1; %number of layers of grading
throwAwayParam = 0; %no need to remove any basis elements
OverSample = 1.25; %choose amount to oversample by (40% here)


% andrew solve
[v_N1, GOA1, colMatrix1, colRHS1, col_points1,...
v_N2, GOA2, colMatrix2, colRHS2, col_points2, VHNA1, VHNA2] ...
    = AG_code_pulling_out_info(pMax, cL, sigmaGrad, nLayers, OverSample, ...
    Gamma1, Gamma2, kwave, uinc );

%% Info needed for our part of the solve
C1 = 1;
C2 = pi;

Lgrad_coeff = 0.15;
alpha = 2;

G1_data.col_points = col_points1;
G2_data.col_points = col_points2;

G1_data.G = [vertices1(1, 1), vertices1(1, 2), vertices1(2, 1), vertices1(2, 2)];
G2_data.G = [vertices2(1, 1), vertices2(1, 2), vertices2(2, 1), vertices2(2, 2)];

G1_data.L = sqrt( (vertices1(2, 1) - vertices1(1, 1))^2 + (vertices1(2, 2) - vertices1(1, 2))^2 );
G1_data.x_col = vertices1(1, 1) + col_points1*(vertices1(2, 1) - vertices1(1, 1))/G1_data.L;
G1_data.y_col = vertices1(1, 2) + col_points1*( vertices1(2, 2) - vertices1(1, 2) )/G1_data.L;

G2_data.L = sqrt( (vertices2(2, 1) - vertices2(1, 1))^2 + (vertices2(2, 2) - vertices2(1, 2))^2 );
G2_data.x_col = vertices2(1, 1) + col_points2*(vertices2(2, 1) - vertices2(1, 1))/G2_data.L;
G2_data.y_col = vertices2(1, 2) + col_points2*( vertices2(2, 2) - vertices2(1, 2) )/G2_data.L;

G1_data.n = [-(G1_data.G(4) - G1_data.G(2)), G1_data.G(3) - G1_data.G(1)]...
    /G1_data.L;
G2_data.n = [-(G2_data.G(4) - G2_data.G(2)), G2_data.G(3) - G2_data.G(1)]...
    /G2_data.L;

G1_data.alpha = -sign(dot(d, G1_data.n));
G2_data.alpha = -sign(dot(d, G2_data.n));
%% Quadrature

G1_data = get_graded_quad_points_HF_it(G1_data, C_wl_quad_outer,...
    C_wl_quad_inner, kwave, Lgrad_coeff, alpha);

G2_data = get_graded_quad_points_HF_it(G2_data, C_wl_quad_outer,...
    C_wl_quad_inner, kwave, Lgrad_coeff, alpha);


%% Iterative solve
[v_N1cell, v_N2cell, phi1_r, phi2_r] = HF_iterative_solve(kwave, ...
    theta, R_max, G1_data, G2_data, VHNA1, colMatrix1, VHNA2, colMatrix2...
    ,vertices1, vertices2, d, C1, C2);
