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
C1 = 1;
C2 = pi;

Lgrad_coeff = 0.15;
alpha = 2;


G1_data.G = [vertices1(1, 1), vertices1(1, 2), vertices1(2, 1), vertices1(2, 2)];
G2_data.G = [vertices2(1, 1), vertices2(1, 2), vertices2(2, 1), vertices2(2, 2)];

G1_data.L = sqrt( (vertices1(2, 1) - vertices1(1, 1))^2 + (vertices1(2, 2) - vertices1(1, 2))^2 );
G1_data.x1_col = vertices1(1, 1) + col_points1*(vertices1(2, 1) - vertices1(1, 1))/G1_data.L;
G1_data.y1_col = vertices1(1, 2) + col_points1*( vertices1(2, 2) - vertices1(1, 2) )/G1_data.L;

G2_data.L = sqrt( (vertices2(2, 1) - vertices2(1, 1))^2 + (vertices2(2, 2) - vertices2(1, 2))^2 );
G2_data.x2_col = vertices2(1, 1) + col_points2*(vertices2(2, 1) - vertices2(1, 1))/G2_data.L;
G2_data.y2_col = vertices2(1, 2) + col_points2*( vertices2(2, 2) - vertices2(1, 2) )/G2_data.L;

G1_data.n = [-(G1_data.G(4) - G1_data.G(2)), G1_data.G(3) - G1_data.G(1)]...
    /G1_data.L;
G2_data.n = [-(G2_data.G(4) - G2_data.G(2)), G2_data.G(3) - G2_data.G(1)]...
    /G2_data.L;

G1_data.alpha = -sign(dot(d, G1_data.n));
G2_data.alpha = -sign(dot(d, G2_data.n));

% quadrature


% outer quadrature, for Sij phij and Sjj Psij
C_wl_quad_outer = 1/20;  % number of quad points per wl
% Inner quadrature, for computing Psij
C_wl_quad_inner = 1/20;  % number of quad points per wl

G1_data = get_graded_quad_points_HF_it(G1_data, C_wl_quad_outer,...
    C_wl_quad_inner, kwave, Lgrad_coeff, alpha);

G2_data = get_graded_quad_points_HF_it(G2_data, C_wl_quad_outer,...
    C_wl_quad_inner, kwave, Lgrad_coeff, alpha);


%% Move onto computing the solve

% 0th iteration, will be able to compare the computation of the RHS to the
% HNABEMLAB code.

% The right hand side in this case is ui - S11duidn evaluated at the
% collocation points
LoB_1_0 = 2*(graded_PIM_int_hankel_f(kwave, col_points1, G1_data.w_outer, ...
    G1_data.t_mid_q_outer, ...
    duidn(vertices1, G1_data.L, kwave, d, G1_data.t_mid_q_outer), ...
    G1_data.t_grid_outer, C1, C2)...
    + graded_PIM_int_hankel_f(kwave, col_points1, flip(G1_data.w_outer), ...
    G1_data.L - flip(G1_data.t_mid_q_outer), ...
    duidn(vertices1, G1_data.L, kwave, d, ...
    G1_data.L - flip(G1_data.t_mid_q_outer)), ...
    G1_data.L - flip(G1_data.t_grid_outer), C1, C2));
RHS1_0 = incident(kwave, theta, G1_data.x1_col, G1_data.y1_col) - LoB_1_0;
%     -(graded_PIM_int_hankel_f(kwave, col_points1, G1_data.w_outer, ...
%     G1_data.t_mid_q_outer, ...
%     duidn(vertices1, G1_data.L, kwave, d, G1_data.t_mid_q_outer), ...
%     G1_data.t_grid_outer, C1, C2));% ...
% 
%     + graded_PIM_int_hankel_f(kwave, col_points1, flip(G1_data.w_outer), ...
%     G1_data.L - flip(G1_data.t_mid_q_outer), ...
%     duidn(vertices1, G1_data.L, kwave, d, ...
%     G1_data.L - flip(G1_data.t_mid_q_outer)), ...
%     G1_data.L - G1_data.t_grid_outer, C1, C2));

% coeff1_0 = 





