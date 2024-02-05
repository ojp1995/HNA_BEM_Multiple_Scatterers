function [aj1, aj2, uinc, G1_data, G2_data] = ...
    compute_2_screen_direct_poly(G1_data, G2_data, k, theta, C1, C2, ...
    C_wl_bf1, C_wl_bf2, C_wl_quad, Lgrad_coeff, alpha)
% In this function we will compute the coefficients for the 2 screen
% problem with the direct polynomial solver


% computing the basis functions
G1_data = get_bf_graded_grid(G1_data, C_wl_bf1, k, Lgrad_coeff, alpha);

G2_data = get_bf_graded_grid(G2_data, C_wl_bf2, k, Lgrad_coeff, ...
    alpha);

% collocation points
G1_data = manipulate_collocation_points_graded(G1_data);
G2_data = manipulate_collocation_points_graded(G2_data);

% quadrature nodes and other info needed
G1_data = get_graded_quad_points(G1_data, C_wl_quad, k, ...
    Lgrad_coeff, alpha);

G2_data = get_graded_quad_points(G2_data, C_wl_quad, k, ...
    Lgrad_coeff, alpha);

% first get the matrices:
[S11, S12, S21, S22, u_inc1, u_inc2] = ...
    compute_matrices_for_iterative_solve(G1_data, G2_data, k, ...
    G1_data.t_bf_grid, G2_data.t_bf_grid, theta, C1, C2 );

% assemble matrix
A = [S11 S12; 
    S21 S22];

uinc = [u_inc1(:); u_inc2(:)];

% solve

aj = A\uinc;

aj1 = aj(1:length(u_inc1));
aj2 = aj(length(u_inc1)+1:end);

