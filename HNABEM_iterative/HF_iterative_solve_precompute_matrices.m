function [v_N1cell, v_N2cell, phi1_r, phi2_r] = ...
    HF_iterative_solve_precompute_matrices(kwave, theta, R_max, ...
    G1_data, G2_data, VHNA1, colMatrix1, VHNA2, colMatrix2...
    ,vertices1, vertices2, d, C1, C2)



%% precompute matrices for iterative solve
[S11, S22, S12, S21, K12, K21, K12_inner, K21_inner] = ...
    compute_matrices_for_HNA_it_solve(G1_data, G2_data, kwave, C1, C2);

% initialising cells for storing the coefficient structures
v_N1cell = {};
v_N2cell = {};
phi1_r = {};
phi2_r = {};

% Step 0 first (special case)
LoB_1_0 = 2*(graded_PIM_int_hankel_f(kwave, G1_data.col_points, ...
    G1_data.w_comb_outer, G1_data.t_mid_q_comb_outer, ...
    duidn(G1_data, G1_data.L, kwave, d, G1_data.t_mid_q_comb_outer), ...
    G1_data.t_grid_comb_outer, C1, C2));

RHS1_0 = incident(kwave, theta, G1_data.x_col, G1_data.y_col) - LoB_1_0;

coeff1_0 = colMatrix1\RHS1_0;

v_N1cell{1} = ProjectionFunction(coeff1_0, VHNA1);

% First compute phi1_0 at outer and inner nodes
v_N1_r = v_N1cell{1};
phi1_r_outer = v_N1_r.eval(G1_data.t_mid_q_comb_outer, 1) ...
    + 2*G1_data.alpha*duidn(G1_data, G1_data.L, kwave, d, ...
    G1_data.t_mid_q_comb_outer);
phi1_r_inner = v_N1_r.eval(G1_data.t_mid_q_comb_inner, 1) ...
    + 2*G1_data.alpha*duidn(G1_data, G1_data.L, kwave, d, ...
    G1_data.t_mid_q_comb_inner);

phi1_r{1} = phi1_r_outer;

S21_phi1_0 = S21*phi1_r_outer;

K21_phi1_0 = K21*(G1_data.beta_inner.*phi1_r_inner);

Psi_2_1 = 2*G2_data.alpha*duidn(G2_data, G2_data.L, kwave, d, ...
    G2_data.t_mid_q_comb_outer) + K21_phi1_0;

S22Psi2_1 = S22*Psi_2_1;

RHS2_1 = incident(kwave, theta, G2_data.x_col, G2_data.y_col) ...
    - S21_phi1_0 - S22Psi2_1;

coeff2_1 = colMatrix2\RHS2_1;

v_N2cell{1}= ProjectionFunction(coeff2_1, VHNA2);

for r = 2:R_max

    % Solving for phi1^(2r-2) iteration
    % Step 1 compute inner and outer phi2^(2r - 3) previous iteration
    % pulling out coefficiens
    v_N2_r = v_N2cell{r - 1};
    phi2_r_outer = v_N2_r.eval(G2_data.t_mid_q_comb_outer, 1) +...
        2*G2_data.alpha*duidn(G2_data, G2_data.L, kwave, d, ...
        G2_data.t_mid_q_comb_outer) ...
        + K21*(G1_data.beta_inner.*phi1_r_inner);

     phi2_r_inner = v_N2_r.eval(G2_data.t_mid_q_comb_inner, 1) +...
        2*G2_data.alpha*duidn(G2_data, G2_data.L, kwave, d, ...
        G2_data.t_mid_q_comb_inner) ...
        + K21_inner*(G1_data.beta_inner.*phi1_r_inner);

      phi2_r{r-1} = phi2_r_outer; % saving for output
    % Compute S12phi2_1, x are the collocation points, integration variables
    % out outer variables
    
    S12_phi2_r = S12*phi2_r_outer;

    K12_phi2_r = K12*(G2_data.beta_inner.*phi2_r_inner);

    Psi_1_r = 2*G1_data.alpha*duidn(G1_data, G1_data.L, kwave, d, ...
        G1_data.t_mid_q_comb_outer) + K12_phi2_r;
    
    S11Psi1_r = S11*Psi_1_r;

    RHS1_r = incident(kwave, theta, G1_data.x_col, G1_data.y_col) ...
        - S12_phi2_r - S11Psi1_r;
    
    % solve
    coeff1_r = colMatrix1\RHS1_r;

    v_N1cell{r} = ProjectionFunction(coeff1_r, VHNA1);

    %% Computing iterative solve for phi2^(2r-1)
    v_N1_r = v_N1cell{r};
    phi1_r_outer = v_N1_r.eval(G1_data.t_mid_q_comb_outer, 1) +...
        2*G1_data.alpha*duidn(G1_data, G1_data.L, kwave, d, ...
        G1_data.t_mid_q_comb_outer) ...
        + K12*(G2_data.beta_inner.*phi2_r_inner);

    phi1_r_inner = v_N1_r.eval(G1_data.t_mid_q_comb_inner, 1) +...
        2*G1_data.alpha*duidn(G1_data, G1_data.L, kwave, d, ...
        G1_data.t_mid_q_comb_inner) ...
        + K12_inner*(G2_data.beta_inner.*phi2_r_inner);

     phi1_r{r} = phi1_r_outer; % saving for output

      % Compute S21phi1_r, x are the collocation points, integration variables
    % out outer variables
    
    S21_phi1_r = S21*phi1_r_outer;

    % Now compute K12 phi2_1, x are G1_outer variables, y are G2_inner
    % variables
    K21_phi1_r = K21*(G1_data.beta_inner.*phi1_r_inner);

    Psi_2_r = 2*G2_data.alpha*duidn(G2_data, G2_data.L, kwave, d, ...
        G2_data.t_mid_q_comb_outer) + K21_phi1_r;

    S22Psi2_r = S22*Psi_2_r;

    RHS2_r = incident(kwave, theta, G2_data.x_col, G2_data.y_col) ...
        - S21_phi1_r - S22Psi2_r;
    
    % solve
    coeff2_r = colMatrix2\RHS2_r;

    v_N2cell{r} = ProjectionFunction(coeff2_r, VHNA2);


end

% computing the final phi2_R
v_N2_r = v_N2cell{r};
phi2_r_outer = v_N2_r.eval(G2_data.t_mid_q_comb_outer, 1) +...
        2*G2_data.alpha*duidn(G2_data, G2_data.L, kwave, d, ...
        G2_data.t_mid_q_comb_outer) ...
        + K21*(G1_data.beta_inner.*phi1_r_inner);


phi2_r{r} = phi2_r_outer; % saving for output











