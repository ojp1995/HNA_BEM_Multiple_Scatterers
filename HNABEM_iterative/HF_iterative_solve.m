function [v_N1cell, v_N2cell, phi1_r, phi2_r] = HF_iterative_solve(kwave, ...
    theta, R_max, G1_data, G2_data, VHNA1, colMatrix1, VHNA2, colMatrix2...
    ,vertices1, vertices2, d, C1, C2)
% In this function we will compute the iterative coefficients and also
% output phi1 and phi2 for a number of iterations
%
% Input parameters

% Output parameters

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

% Step 1 - Doing outside for loop for ease (and as an example for the 
% process)
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
v_N1_r = v_N1cell{1};
phi1_0_outer = v_N1_r.eval(G1_data.t_mid_q_comb_outer, 1) ...
    + 2*G1_data.alpha*duidn(G1_data, G1_data.L, kwave, d, G1_data.t_mid_q_comb_outer);
phi1_0_inner = v_N1_r.eval(G1_data.t_mid_q_comb_inner, 1) ...
    + 2*G1_data.alpha*duidn(G1_data, G1_data.L, kwave, d, G1_data.t_mid_q_comb_inner);

phi1_r{1} = phi1_0_outer;
% now lets compute S21 phi1_0
S21_phi1_0 = midpoint_hankel_f_diff_screen(kwave, G2_data.x_col, ...
    G2_data.y_col, G1_data.x_q_comb_outer, G1_data.y_q_comb_outer, ...
    G1_data.w_comb_outer, phi1_0_outer);


% Now compute Psi2_1 evaluated at G2_data.t_mid_q_comb_outer using the
% inner nodes for integration

K21_phi1_0 =  G1_data.beta_outer.*midpoint_dphikdn_f_diff_screen(kwave, ...
    G2_data.x_q_comb_outer, G2_data.y_q_comb_outer, ...
    G1_data.w_comb_inner, G1_data.x_q_comb_inner, ...
    G1_data.y_q_comb_inner, phi1_0_inner, G2_data.n);

Psi_2_1 = 2*G2_data.alpha*duidn(G2_data, G2_data.L, kwave, d, ...
    G2_data.t_mid_q_comb_outer) + K21_phi1_0;

% Now compute S22Psi2_1
S22Psi2_1 = graded_PIM_int_hankel_f(kwave, G2_data.col_points,...
    G2_data.w_comb_outer, G2_data.t_mid_q_comb_outer, Psi_2_1, ...
    G2_data.t_grid_comb_outer, C1, C2);

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
        + G1_data.beta_outer.*midpoint_dphikdn_f_diff_screen(kwave, ...
        G2_data.x_q_comb_outer, G2_data.y_q_comb_outer, ...
        G1_data.w_comb_inner, G1_data.x_q_comb_inner, ...
        G1_data.y_q_comb_inner, phi1_0_inner, G2_data.n);

    phi2_r_inner = v_N2_r.eval(G2_data.t_mid_q_comb_inner, 1) +...
        2*G2_data.alpha*duidn(G2_data, G2_data.L, kwave, d, ...
        G2_data.t_mid_q_comb_inner) ...
        + G1_data.beta_inner.*midpoint_dphikdn_f_diff_screen(kwave, ...
        G2_data.x_q_comb_inner, G2_data.y_q_comb_inner, ...
        G1_data.w_comb_inner, G1_data.x_q_comb_inner, ...
        G1_data.y_q_comb_inner, phi1_0_inner, G2_data.n);
    
    phi2_r{r-1} = phi2_r_outer; % saving for output
    % Compute S12phi2_1, x are the collocation points, integration variables
    % out outer variables
    
    S12_phi2_r = midpoint_hankel_f_diff_screen(kwave, G1_data.x_col, ...
        G1_data.y_col, G2_data.x_q_comb_outer, G2_data.y_q_comb_outer, ...
        G2_data.w_comb_outer, phi2_r_outer);
    % Now compute K12 phi2_1, x are G1_outer variables, y are G2_inner
    % variables
    K12_phi2_r = G2_data.beta_outer.*midpoint_dphikdn_f_diff_screen(kwave, ...
        G1_data.x_q_comb_outer, G1_data.y_q_comb_outer, ...
        G2_data.w_comb_inner, G2_data.x_q_comb_inner, ...
        G2_data.y_q_comb_inner, phi2_r_inner, G1_data.n);
    
    
    Psi_1_r = 2*G1_data.alpha*duidn(G1_data, G1_data.L, kwave, d, ...
        G1_data.t_mid_q_comb_outer) + K12_phi2_r;
    
    S11Psi1_r = graded_PIM_int_hankel_f(kwave, G1_data.col_points,...
        G1_data.w_comb_outer, G1_data.t_mid_q_comb_outer, Psi_1_r, ...
        G1_data.t_grid_comb_outer, C1, C2);
    
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
        + G2_data.beta_outer.*midpoint_dphikdn_f_diff_screen(kwave, ...
        G1_data.x_q_comb_outer, G1_data.y_q_comb_outer, ...
        G2_data.w_comb_inner, G2_data.x_q_comb_inner, ...
        G2_data.y_q_comb_inner, phi2_r_inner, G1_data.n);

    phi1_r_inner = v_N1_r.eval(G1_data.t_mid_q_comb_inner, 1) +...
        2*G1_data.alpha*duidn(G1_data, G1_data.L, kwave, d, ...
        G1_data.t_mid_q_comb_inner) ...
        + G2_data.beta_inner.*midpoint_dphikdn_f_diff_screen(kwave, ...
        G1_data.x_q_comb_inner, G1_data.y_q_comb_inner, ...
        G2_data.w_comb_inner, G2_data.x_q_comb_inner, ...
        G2_data.y_q_comb_inner, phi2_r_inner, G1_data.n);
    
    phi1_r{r} = phi1_r_outer; % saving for output

   
    % Compute S21phi1_r, x are the collocation points, integration variables
    % out outer variables
    
    S21_phi1_r = midpoint_hankel_f_diff_screen(kwave, G2_data.x_col, ...
        G2_data.y_col, G1_data.x_q_comb_outer, G1_data.y_q_comb_outer, ...
        G1_data.w_comb_outer, phi1_r_outer);
    
    % Now compute K12 phi2_1, x are G1_outer variables, y are G2_inner
    % variables
    K21_phi1_r =  G1_data.beta_outer.*midpoint_dphikdn_f_diff_screen(kwave, ...
        G2_data.x_q_comb_outer, G2_data.y_q_comb_outer, ...
        G1_data.w_comb_inner, G1_data.x_q_comb_inner, ...
        G1_data.y_q_comb_inner, phi1_r_inner, G2_data.n);

    Psi_2_r = 2*G2_data.alpha*duidn(G2_data, G2_data.L, kwave, d, ...
        G2_data.t_mid_q_comb_outer) + K21_phi1_r;
    
    S22Psi2_r = graded_PIM_int_hankel_f(kwave, G2_data.col_points,...
        G2_data.w_comb_outer, G2_data.t_mid_q_comb_outer, Psi_2_r, ...
        G2_data.t_grid_comb_outer, C1, C2);
    
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
        + G1_data.beta_outer.*midpoint_dphikdn_f_diff_screen(kwave, ...
        G2_data.x_q_comb_outer, G2_data.y_q_comb_outer, ...
        G1_data.w_comb_inner, G1_data.x_q_comb_inner, ...
        G1_data.y_q_comb_inner, phi1_0_inner, G2_data.n);

    
    phi2_r{r} = phi2_r_outer; % saving for output



