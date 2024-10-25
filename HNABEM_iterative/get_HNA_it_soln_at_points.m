function [phi1, phi2] = get_HNA_it_soln_at_points(G1_data, G2_data, k, ...
    d, v1_coeffs, v2_coeffs, Rmax)
% Evaluating at points phi1_points and phi2_points


phi1 = zeros(length(G1_data.t_mid_q_comb_outer), length(Rmax));
phi2 = zeros(length(G2_data.t_mid_q_comb_outer), length(Rmax));

v_N1_r = v1_coeffs{1};

phi1(:, 1) = v_N1_r.eval(G1_data.t_mid_q_comb_outer, 1) ...
    + 2*G1_data.alpha*duidn(G1_data, G1_data.L, k, d, ...
    G1_data.t_mid_q_comb_outer);

phi1_HNA_eval_inner = zeros(length(G1_data.t_mid_q_comb_inner), Rmax);
phi2_HNA_eval_inner = zeros(length(G2_data.t_mid_q_comb_inner), Rmax);

phi1_HNA_eval_inner(:, 1) = v_N1_r.eval(G1_data.t_mid_q_comb_inner, 1) ...
    + 2*G1_data.alpha*duidn(G1_data, G1_data.L, k, d, ...
    G1_data.t_mid_q_comb_inner);

for r = 2:Rmax

    v_N2_r = v2_coeffs{r - 1};

    % compute phi2 outer
    phi2(:, r - 1) = v_N2_r.eval( G2_data.t_mid_q_comb_outer, 1) +...
        2*G2_data.alpha*duidn(G2_data, G2_data.L, k, d, ...
        G2_data.t_mid_q_comb_outer) ...
        + midpoint_dphikdn_f_diff_screen(k, ...
        G2_data.x_q_comb_outer, G2_data.y_q_comb_outer, ...
        G1_data.w_comb_inner, G1_data.x_q_comb_inner, ...
        G1_data.y_q_comb_inner, G1_data.beta_inner.*phi1_HNA_eval_inner(:, r - 1), ...
        G2_data.n);

    % compute phi2 inner
    phi2_HNA_eval_inner(:, r - 1) = v_N2_r.eval(G2_data.t_mid_q_comb_inner, 1) +...
        2*G2_data.alpha*duidn(G2_data, G2_data.L, k, d, ...
        G2_data.t_mid_q_comb_inner) ...
        + midpoint_dphikdn_f_diff_screen(k, ...
        G2_data.x_q_comb_inner, G2_data.y_q_comb_inner, ...
        G1_data.w_comb_inner, G1_data.x_q_comb_inner, ...
        G1_data.y_q_comb_inner, G1_data.beta_inner.*phi1_HNA_eval_inner(:, r-1), ...
        G2_data.n);

    % computing the phi1 outer

    v_N1_r = v1_coeffs{r};


    phi1(:, r) = v_N1_r.eval(G1_data.t_mid_q_comb_outer, 1) +...
        2*G1_data.alpha*duidn(G1_data, G1_data.L, k, d, ...
        G1_data.t_mid_q_comb_outer) ...
        + midpoint_dphikdn_f_diff_screen(k, ...
        G1_data.x_q_comb_outer, G1_data.y_q_comb_outer, ...
        G2_data.w_comb_inner, G2_data.x_q_comb_inner, ...
        G2_data.y_q_comb_inner, G2_data.beta_inner.*phi2_HNA_eval_inner(:, r - 1), ...
        G1_data.n);

 % phi1 inner
     phi1_HNA_eval_inner(:, r) = v_N1_r.eval(G1_data.t_mid_q_comb_inner, 1) +...
        2*G1_data.alpha*duidn(G1_data, G1_data.L, k, d, ...
        G1_data.t_mid_q_comb_inner) ...
        + midpoint_dphikdn_f_diff_screen(k, ...
        G1_data.x_q_comb_inner, G1_data.y_q_comb_inner, ...
        G2_data.w_comb_inner, G2_data.x_q_comb_inner, ...
        G2_data.y_q_comb_inner, G2_data.beta_inner.*phi2_HNA_eval_inner(:, r - 1), ...
        G1_data.n);



    

end

% final phi2 calc
v_N2_r = v2_coeffs{end};

phi2(:, end) = v_N2_r.eval(G2_data.t_mid_q_comb_outer, 1) +...
        2*G2_data.alpha*duidn(G2_data, G2_data.L, k, d, ...
        G2_data.t_mid_q_comb_outer) ...
        + midpoint_dphikdn_f_diff_screen(k, ...
       G2_data.x_q_comb_outer, G2_data.y_q_comb_outer, ...
        G1_data.w_comb_inner, G1_data.x_q_comb_inner, ...
        G1_data.y_q_comb_inner, G1_data.beta_inner.*phi1_HNA_eval_inner(:, end), ...
        G2_data.n);