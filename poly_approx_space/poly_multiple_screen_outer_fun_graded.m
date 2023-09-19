% outer function for the polynomial solver with graded quadrature
clear all

% only adding one path to centralise solvers
addpath('../General_functions/')

% switch and things get interesting
G1_data.G = [-2*pi, 2*pi, 0, 0];
% G1_data.L = sqrt( (G1(3) - G1(1))^2 +(G1(4) - G1(2))^2 );

G2_data.G = [2*pi, 0, 5*pi, 3*pi]; 
% G2_data.L = sqrt( (G2(3) - G2(1))^2 +(G2(4) - G2(2))^2 );

Lgrad_coeff = 0.15;
alpha = 4;

C_wl= 1/4;

k = 1;  % wavenumber

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
% h1 = L1/1000;
% h2 = L2/1000;
% t1_grid = [0: h1: L1];
% t2_grid = [0: h2: L2];

% The support for the basis functions

% C_wl_bf1 = 1/20;
% C_wl_bf2 = 1/20;
% [~, ~, ~, ~, t1_bf_grid, ~, ~, ~, ~] = discretistion_vars_graded(...
%     G1_data.G, C_wl_bf1, k, Lgrad_coeff, alpha);
% [~, ~, ~, ~, t2_bf_grid, ~, ~, ~, ~] = discretistion_vars_graded(...
%     G2_data.G, C_wl_bf2, k, Lgrad_coeff, alpha);

% quadrature nodes and other information needed
[G1_data.x_1_q, G1_data.y_1_q, G1_data.x_2_q, G1_data.y_2_q, G1_data.t_grid, G1_data.t_mid, G1_data.w, G1_data.N, G1_data.L] = ...
    discretistion_vars_graded(G1_data.G, C_wl, k, Lgrad_coeff, alpha);

[G2_data.x_1_q, G2_data.y_1_q, G2_data.x_2_q, G2_data.y_2_q, G2_data.t_grid, G2_data.t_mid, G2_data.w, G2_data.N, G2_data.L] = ...
    discretistion_vars_graded(G2_data.G, C_wl, k, Lgrad_coeff, alpha);


t1_bf_grid = linspace(0, G1_data.L/2, 10);
t2_bf_grid = linspace(0, G2_data.L/2, 10);

% Collocation points (using the same grid currently but will change)
% col_choice1 = sort(randi(length(G1_data.t_mid(:)), 20, 1));

G1_data.s = [ G1_data.t_mid(1:end) ; flip(G1_data.L - G1_data.t_mid(1:end)) ];
G1_data.x_col = [ G1_data.x_1_q(1:end) ; G1_data.x_2_q(1:end) ];
G1_data.y_col = [ G1_data.y_1_q(1:end) ; G1_data.y_2_q(1:end) ];

% col_choice2 = sort(randi(length(G2_data.t_mid(:)), 20, 1));
G2_data.s = [ G2_data.t_mid(1:end) ; flip(G2_data.L - G2_data.t_mid(1:end)) ];
G2_data.x_col = [ G2_data.x_1_q(1:end) ; G2_data.x_2_q(1:end) ];
G2_data.y_col = [ G2_data.y_1_q(1:end) ; G2_data.y_2_q(1:end) ];

[S11, S12, S21, S22, u_inc1, u_inc2] = ...
    compute_matrices_for_iterative_solve(G1_data, G2_data, k, ...
    t1_bf_grid, t2_bf_grid, theta, C1, C2 );

% all in one solve
A = [S11  S12 ; S21 S22];

u_inc = [u_inc1 ; u_inc2];

coeffs = A\u_inc;

aj_1 = coeffs(1:2*length(t1_bf_grid) - 2);
aj_2 = coeffs(2*length(t1_bf_grid) - 2+ 1: end);

%%
% now we need to find a way to plot these coefficients and also make sure
% we don't lose information if they are too close to the edge
x1_plot = linspace(0.01, G1_data.L/2 - 0.01, 200);

x1_plotting = [x1_plot(:) ; (G1_data.L/2 + x1_plot(:) )];

phi1 = graded_coeff_2_solution(S11\u_inc1, t1_bf_grid, x1_plot, G1_data.L);

figure()
plot(x1_plotting, phi1)