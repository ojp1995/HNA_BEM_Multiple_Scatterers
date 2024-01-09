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

% quadrature


% outer quadrature, for Sij phij and Sjj Psij
C_wl_quad_outer = 1/20;  % number of quad points per wl
% Inner quadrature, for computing Psij
C_wl_quad_inner = 1/30;  % number of quad points per wl

G1_data = get_graded_quad_points_HF_it(G1_data, C_wl_quad_outer,...
    C_wl_quad_inner, kwave, Lgrad_coeff, alpha);

G2_data = get_graded_quad_points_HF_it(G2_data, C_wl_quad_outer,...
    C_wl_quad_inner, kwave, Lgrad_coeff, alpha);


%% Move onto computing the solve

% 0th iteration, will be able to compare the computation of the RHS to the
% HNABEMLAB code.

% The right hand side in this case is ui - S11duidn evaluated at the
% collocation points
LoB_1_0 = 2*(graded_PIM_int_hankel_f(kwave, col_points1, ...
    G1_data.w_comb_outer, G1_data.t_mid_q_comb_outer, ...
    duidn(vertices1, G1_data.L, kwave, d, G1_data.t_mid_q_comb_outer), ...
    G1_data.t_grid_comb_outer, C1, C2));

% % % % LoB_1_0_test = 2*( graded_PIM_int_hankel_f(kwave, col_points1, G1_data.w_outer, ...
% % % %     G1_data.t_mid_q_outer, ...
% % % %     duidn(vertices1, G1_data.L, kwave, d, G1_data.t_mid_q_outer), ...
% % % %     G1_data.t_grid_outer, C1, C2)...
% % % %     + graded_rescalled_PIM_int_hankel_f(kwave, col_points1,...
% % % %     G1_data.w_outer, ...
% % % %     G1_data.t_mid_q_outer, ...
% % % %     duidn(vertices1, G1_data.L, kwave, d, G1_data.t_mid_q_outer), ...
% % % %     G1_data.t_grid_outer, C1, C2, G1_data.L));
% 2*(graded_PIM_int_hankel_f(kwave, col_points1, G1_data.w_outer, ...
%     G1_data.t_mid_q_outer, ...
%     duidn(vertices1, G1_data.L, kwave, d, G1_data.t_mid_q_outer), ...
%     G1_data.t_grid_outer, C1, C2)...
%     + graded_PIM_int_hankel_f(kwave, col_points1, flip(G1_data.w_outer), ...
%     G1_data.L - flip(G1_data.t_mid_q_outer), ...
%     duidn(vertices1, G1_data.L, kwave, d, ...
%     G1_data.L - flip(G1_data.t_mid_q_outer)), ...
%     G1_data.L - flip(G1_data.t_grid_outer), C1, C2));
RHS1_0 = incident(kwave, theta, G1_data.x_col, G1_data.y_col) - LoB_1_0;
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

% quick test, Mostly the same, there are some differences, could be solved
% with convergence/ rescaling integrals I believe
% [colRHS1, RHS1_0]  

%% solve
coeff1_0 = colMatrix1\RHS1_0;

v_N1_0 = ProjectionFunction(coeff1_0, VHNA1);

[coeff1_0, v_N1_0.coeffs(), v_N1.coeffs()]

%% construct the solution
phi1_0 = @(x) v_N1_0.eval(x, 1) + 2*duidn(vertices1, G1_data.L, kwave, ...
    d, x);

figure();
plot(G1_data.t_mid_q_comb_outer/G1_data.L, phi1_0(G1_data.t_mid_q_comb_outer))
hold on
plot(G1_data.t_mid_q_comb_outer/G1_data.L, ...
    v_N1.eval(G1_data.t_mid_q_comb_outer, 1) ...
    + GOA1.eval(G1_data.t_mid_q_comb_outer, 1), '--')
xlim([-0.05 1.05])
ylim([-30 30])


%% Step 1
% We have three parts to compute:
    % ui (x, col points)
    % S21phi1_0 (x, col points. phi1_0 evaluated at the outer nodes (y))
    % S22 Psi2_1  (x, col poionts, Psi2_1 evaluated at the outer nodes (y))
    %
    % where Psi2_1 can be broken down as: (x here is outer quadrature points)
        % 2duidn (quadrature points from outer integral)
        % S21 phi1_0 (x, quadrature points from outer integral, 
                        % phi needs to be evaluated at inner nodes (y))

% First compute phi1_0 at outer and inner nodes
phi1_0_outer = phi1_0(G1_data.t_mid_q_comb_outer);
phi1_0_inner = phi1_0(G1_data.t_mid_q_comb_inner);



% now lets compute S21 phi1_0
S21_phi1_0 = midpoint_hankel_f_diff_screen(kwave, G2_data.x_col, ...
    G2_data.y_col, G1_data.x_q_comb_outer, G1_data.y_q_comb_outer, ...
    G1_data.w_comb_outer, phi1_0_outer);

figure(); plot(col_points2/G2_data.L, real(S21_phi1_0)); title('S21 $\phi_1^0$ evaluated at collocation points')
% debugging
S21_phi1_0_test = midpoint_hankel_f_diff_screen(kwave, G2_data.x_q_comb_outer, ...
G2_data.y_q_comb_outer, G1_data.x_q_comb_outer, G1_data.y_q_comb_outer, ...
 G1_data.w_comb_outer, phi1_0_outer);

figure(); plot(G2_data.t_mid_q_comb_outer /G2_data.L, real(S21_phi1_0_test)); title('Hybrid S21 $\phi_1^0$ evaluated at more points')




% Now compute Psi2_1 evaluated at G2_data.t_mid_q_comb_outer using the
% inner nodes for integration

K21_phi1_0 =  midpoint_dphikdn_f_diff_screen(kwave, ...
    G2_data.x_q_comb_outer, G2_data.y_q_comb_outer, ...
    G1_data.w_comb_inner, G1_data.x_q_comb_inner, ...
    G1_data.y_q_comb_inner, phi1_0_inner, G2_data.n);

quad_mid = length(G1_data.w_comb_inner)/2;
K21_phi1_0_half1_nonscaled =  midpoint_dphikdn_f_diff_screen(kwave, ...
    G2_data.x_q_comb_outer, G2_data.y_q_comb_outer, ...
    G1_data.w_comb_inner(1:quad_mid), G1_data.x_q_comb_inner(1:quad_mid), ...
    G1_data.y_q_comb_inner(1:quad_mid), phi1_0_inner(1:quad_mid), G2_data.n);

K21_phi1_0_half2_nonscaled =  midpoint_dphikdn_f_diff_screen(kwave, ...
    G2_data.x_q_comb_outer, G2_data.y_q_comb_outer, ...
    G1_data.w_comb_inner(quad_mid+1:end), G1_data.x_q_comb_inner(quad_mid+1:end), ...
    G1_data.y_q_comb_inner(quad_mid+1:end), phi1_0_inner(quad_mid+1:end), G2_data.n);

K21_phi1_0_half1 = midpoint_dphikdn_f_diff_screen(kwave,...
    G2_data.x_q_comb_outer, G2_data.y_q_comb_outer, ...
    G1_data.w_inner, G1_data.x_1_q_inner, ...
    G1_data.y_1_q_inner, ...
    phi1_0(G1_data.t_mid_q_inner), G2_data.n);

K21_phi1_0_half2 = midpoint_dphikdn_f_diff_screen(kwave,...
    G2_data.x_q_comb_outer, G2_data.y_q_comb_outer, ...
    G1_data.w_inner, G1_data.L - G1_data.x_1_q_inner, ...
    G1_data.L - G1_data.y_1_q_inner, ...
    phi1_0(G1_data.L - G1_data.t_mid_q_inner), G2_data.n);

K21phi1_0_half2_mannually_scalled = ...
    recalled_midpoint_dphikdn_f_diff_screen(kwave, ...
    G2_data.x_q_comb_outer, G2_data.y_q_comb_outer, ...
    G1_data.w_inner, G1_data.x_1_q_inner, ...
    G1_data.y_1_q_inner, ...
    phi1_0(G1_data.L - G1_data.t_mid_q_inner), G2_data.n, G1_data.L);

%% Not very convinced lets see how matlab solves it
% anon functions for integration variables and kernel
G1x = @(t) G1_data.G(1) + t.*(G1_data.G(3) - G1_data.G(1))/G1_data.L;
G1y = @(t) G1_data.G(2) + t.*(G1_data.G(4) - G1_data.G(2))/G1_data.L;

dphik_kernel = @(t) 1i*kwave*besselh(1, 1, kwave*sqrt( (G2_data.x_q_comb_outer  - ...
    G1x(t)).^2 + (G2_data.x_q_comb_outer - G1y(t)).^2 ) ).*((G2_data.x_q_comb_outer  - ...
    G1x(t))*G2_data.n(1) + (G2_data.x_q_comb_outer - G1y(t))*G2_data.n(2))...
    .*phi1_0(t)./(2* ...
    sqrt( (G2_data.x_q_comb_outer  -  G1x(t)).^2 + (G2_data.x_q_comb_outer - G1y(t)).^2) );


mat_K21phi1_0 = integral(@(t) dphik_kernel(t), 0, G1_data.L, 'ArrayValued', true);
%% comparison
figure(); 
plot(G2_data.t_mid_q_comb_outer/G2_data.L, real(K21_phi1_0_half1_nonscaled), 'DisplayName', 'Non-scaled, 1st half');
hold on
plot(G2_data.t_mid_q_comb_outer/G2_data.L, real(K21_phi1_0_half1), '--', 'DisplayName', 'scaled, 1st half');
plot(G2_data.t_mid_q_comb_outer/G2_data.L, real(mat_K21phi1_0), '-.', 'DisplayName', 'Matlab 1st half interval')
figure(); 
plot(G2_data.t_mid_q_comb_outer/G2_data.L, real(K21_phi1_0_half2_nonscaled), 'DisplayName', 'Non-scaled, 2nd half');
hold on
plot(G2_data.t_mid_q_comb_outer/G2_data.L, real(K21_phi1_0_half2), '--', 'DisplayName', 'scaled, 2nd half half');
plot(G2_data.t_mid_q_comb_outer/G2_data.L, K21phi1_0_half2_mannually_scalled, '-.', 'DisplayName', 'manually scalled')
plot(G2_data.t_mid_q_comb_outer/G2_data.L, real(mat_K21phi1_0), '-.', 'DisplayName', 'Matlab 2nd half interval')
legend show
Psi_2_1 = 2*G2_data.alpha*duidn(vertices2, G2_data.L, kwave, d, ...
    G2_data.t_mid_q_comb_outer) + K21_phi1_0; %+ K21_phi1_0_half1 + K21_phi1_0_half2;

% Now compute S22Psi2_1
S22Psi2_1 = graded_PIM_int_hankel_f(kwave, col_points2,...
    G2_data.w_comb_outer, G2_data.t_mid_q_comb_outer, Psi_2_1, ...
    G2_data.t_grid_comb_outer, C1, C2);

RHS2_1 = incident(kwave, theta, G2_data.x_col, G2_data.y_col) ...
    - S21_phi1_0 - S22Psi2_1;
figure();
plot(col_points2/G2_data.L, RHS2_1 )

%% solve
coeff2_1 = colMatrix2\RHS2_1;

v_N2_1 = ProjectionFunction(coeff2_1, VHNA2);

%% compute solution
% phi2_1 = @(x) v_N2_1.eval(x, 1) 

phi2_1 = v_N2_1.eval(G2_data.t_mid_q_comb_outer, 1) +...
    2*G2_data.alpha*duidn(vertices2, G2_data.L, kwave, d, ...
    G2_data.t_mid_q_comb_outer) + K21_phi1_0;

figure()
plot(G2_data.t_mid_q_comb_outer/G2_data.L, real(phi2_1))
xlim([-0.05 1.05])
ylim([-30 30])











