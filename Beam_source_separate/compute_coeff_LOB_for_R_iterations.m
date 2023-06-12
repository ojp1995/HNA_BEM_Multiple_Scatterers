function [aj_1_r, aj_2_r, phi1_r_outer, phi2_r_outer] = ...
    compute_coeff_LOB_for_R_iterations(kwave, N_approx, G1, G2,...
    vertices1, vertices2, R, theta,...
    col_points1, x1_col, y1_col, col_points2, x2_col, y2_col, VHNA1,...
    VHNA2, colMatrix1, colMatrix2, d, n1, n2, C1, C2)
% In this function we will be computing the solution for R iterations for
% our two screen problem
%
%
% Input parameters:
%
% Output parameters:
% 

% First we want to compute the numerical parameters:
[x1, y1, t1, t1_mid, h1, h1vector, N1, L1] = ...
    discretisation_variables(G1, N_approx, kwave);
[x2, y2, t2, t2_mid, h2, h2vector, N2, L2] = ...
    discretisation_variables(G2, N_approx, kwave);

% variables for the inner intergals
N_approx_inner = N_approx/2;
[y1nq_1_inner, y2nq_1_inner, ~, t1_mid_inner, h1_inner, ~, ~, ~] = ...
    discretisation_variables(G1, N_approx_inner, kwave);
[y1nq_2_inner, y2nq_2_inner, ~, t2_mid_inner, h2_inner, ~, ~, ~] = ...
    discretisation_variables(G2, N_approx_inner, kwave);

phi1_r_outer = zeros(R+1, length(t1_mid));
phi1_r_inner = zeros(R+1, length(t1_mid_inner));


phi2_r_outer = zeros(R, length(t2_mid));
phi2_r_inner = zeros(R, length(t2_mid_inner));

% Now computing step 0

[v_N_G1_r0] = multiple_Scattering_2screen_step0(kwave, theta, ...
    d, vertices1, L1, col_points1, x1_col, y1_col, h1, t1_mid, t1, colMatrix1, VHNA1, C1, C2);

aj_1_r = {v_N_G1_r0}; % storing the coefficient object in a cell array

phi1_0 = @(x) v_N_G1_r0.eval(x, 1) + 2*duidn(vertices1, L1, kwave, d, x);

phi1_r_outer(1, :) = phi1_0(t1_mid.');
phi1_r_inner(1, :) = phi1_0(t1_mid_inner.');
% compute step 1 (this may get absorbed into the below for loop)
[f_2_2rp1(1, :), f_2_2rp1_uinc(1, :), f_2_2rp1_beam(1, :), ~,~, ~, ~,~, ~] ...
    = beam_compute_RHS_vec_given_coll_vec(vertices2, L2, kwave, d, ...
        theta, n1, col_points2, x2_col, y2_col, C1, C2, x1, y1, ...
        h1, phi1_r_outer(1, :).', x2, y2, t2_mid, ...
        t2, h2, y1nq_1_inner, y2nq_1_inner, h1_inner, ...
        phi1_r_inner(1, :).');
    
v_N_G2_r1 = compute_coeffs_given_A_and_f(colMatrix2, f_2_2rp1(1, :).', VHNA2);

aj_2_r = {v_N_G2_r1};

phi2_r_outer(1, :) = get_phi_j_r(v_N_G2_r1, vertices2, L2, kwave, d, h1, x1, ...
y1, n2, t2_mid, x2, y2, phi1_r_outer(1, :).');

phi2_r_inner(1, :) = get_phi_j_r(v_N_G2_r1, vertices2, L2, kwave, d, h1, x1, ...
y1, n2, t2_mid_inner, y1nq_2_inner, y2nq_2_inner, phi1_r_outer(1, :).');


for r = 1:R
    
    % computing \phi1_2r
    [f_1_2r(r+1,:), f_1_2_uinc(r+1, :), f_1_2_beam(r+1, :),~,~,~,~,~,~] ...
    = beam_compute_RHS_vec_given_coll_vec(vertices1, L1, kwave, d, ...
        theta, n2, col_points1, x1_col, y1_col, C1, C2, x2, y2, ...
        h2, phi2_r_outer(r, :).', x1, y1, t1_mid, ...
        t1, h1, y1nq_2_inner, y2nq_2_inner, h2_inner, ...
        phi2_r_inner(r, :).');
    
    v_N_G1_2r = ...
        compute_coeffs_given_A_and_f(colMatrix1, f_1_2r(r+1, :).', VHNA1);
    
    aj_1_r{r+1} = v_N_G1_2r;
    
    phi1_r_outer(r+1, :) = get_phi_j_r(v_N_G1_2r, vertices1, L1, kwave, d, h2, x2, ...
    y2, n1, t1_mid, x1, y1, phi2_r_outer(r, :).');


    phi1_r_inner(r+1, :) = get_phi_j_r(v_N_G1_2r, vertices1, L1, kwave, d, h2, x2, ...
    y2, n1, t1_mid_inner, y1nq_1_inner, y2nq_1_inner, phi2_r_outer(r, :).');


    % computing \phi2_2r+1
    [f_2_2rp1(r+1, :), f_2_2rp1_uinc(r+1, :), f_2_2rp1_beam(r+1, :), ~,~, ~, ~,~, ~] ...
        = beam_compute_RHS_vec_given_coll_vec(vertices2, L2, kwave, d, ...
        theta, n1, col_points2, x2_col, y2_col, C1, C2, x1, y1, ...
        h1, phi1_r_outer(r+1, :).', x2, y2, t2_mid, ...
        t2, h2, y1nq_1_inner, y2nq_1_inner, h1_inner, ...
         phi1_r_inner(r+1, :).');
     
    v_N_G2_2rp1 = compute_coeffs_given_A_and_f(colMatrix2, f_2_2rp1(r+1, :).', VHNA2);
     
    aj_2_r{r+1} = v_N_G2_2rp1;
     
    phi2_r_outer(r+1, :) = get_phi_j_r(v_N_G2_r1, vertices2, L2, kwave, d, h1, x1, ...
    y1, n2, t2_mid, x2, y2, phi1_r_outer(1, :).');
    
    phi2_r_inner(r+1, :) = get_phi_j_r(v_N_G2_r1, vertices2, L2, kwave, d, h1, x1, ...
    y1, n2, t2_mid_inner, y1nq_2_inner, y2nq_2_inner, phi1_r_outer(1, :).');

     


end


end

