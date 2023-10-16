function [aj_1R, aj_2R, phi_1r, phi_2r] = ...
    iterative_poly_graded_PIM_solve(S11, S12, S21, S22, u_inc1, u_inc2, ...
    R_max, t1_bf_grid, t2_bf_grid, G1_data, G2_data)
% In this function we will be computing the iterative solve for the
% polynomial solver using the graded polynomial graded PIM and graded
% midpoint rule
%
%

aj_1R = zeros(size(S11, 1), 1);
aj_2R = zeros(size(S22, 1), 1);

phi_1r = zeros(length(G1_data.s), 1);
phi_2r = zeros(length(G2_data.s), 1);

aj_1R(:, 1) = S11\u_inc1;

phi_1r(:, 1) = graded_coeff_2_solution(aj_1R(:, 1), t1_bf_grid, ...
    G1_data.t_mid_col, G1_data.L);

aj_2R(:, 1) = S22\(u_inc2 - S21*phi_1r(:, 1));

phi_2r(:, 1) = graded_coeff_2_solution(aj_2R(:, 1), t2_bf_grid, ...
    G2_data.t_mid_col, G2_data.L);
for r = 2:R_max

    aj_1R(:, r) = S11\(u_inc1 - S12*phi_2r(:, r-1));
    
    phi_1r(:, r) = graded_coeff_2_solution(aj_1R(:, r), t1_bf_grid, ...
    G1_data.t_mid_col, G1_data.L);

    aj_2R(:, r) = S22\(u_inc2 - S21*phi_1r(:, r));

    phi_2r(:, r) = graded_coeff_2_solution(aj_2R(:, r), t2_bf_grid, ...
        G2_data.t_mid_col, G2_data.L);


end