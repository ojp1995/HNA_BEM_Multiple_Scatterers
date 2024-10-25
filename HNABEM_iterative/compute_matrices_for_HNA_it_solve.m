function [S11, S22, S12, S21, K12, K21, K12_inner, K21_inner] = ...
    compute_matrices_for_HNA_it_solve(G1_data, G2_data, k, C1, C2)
% This funciton will compute the matricies up front for operators ready to
% then be multiplied by the relevant vectors. S11 and S22, S12 and S21 
% are evaluated at the collocation points. K12, K21,  are for the inner 
% points

% outer matrix S12
S12 = zeros(length(G1_data.col_points), length(G2_data.t_mid_q_comb_outer));
S11 = zeros(length(G1_data.col_points), length(G1_data.t_mid_q_comb_outer));
for j = 1:length(G1_data.col_points)

    dist = sqrt( (G1_data.x_col(j) - G2_data.x_q_comb_outer).^2 ...
        +  (G1_data.y_col(j) - G2_data.y_q_comb_outer).^2 );

    S12(j, :) = 1i.*G2_data.w_comb_outer.*besselh(0, k*dist(:))./4;

    % S11 computation
    t_lower = G1_data.t_grid_comb_outer(1:end - 1);
    t_upper = G1_data.t_grid_comb_outer(2:end);

    S11(j, :) = m1(k, G1_data.col_points(j), G1_data.t_mid_q_comb_outer, C1, C2)...
        .*w1_weights(k, G1_data.col_points(j), t_lower, t_upper) ...
        + G1_data.w_comb_outer.*m2(k, G1_data.col_points(j), ...
        G1_data.t_mid_q_comb_outer, C1, C2 );


end

% outer matrix S21
S21 = zeros(length(G2_data.col_points), length(G1_data.t_mid_q_comb_outer));
S22 = zeros(length(G2_data.col_points), length(G2_data.t_mid_q_comb_outer));
for j = 1:length(G2_data.col_points)

    dist = sqrt( (G2_data.x_col(j) - G1_data.x_q_comb_outer).^2 ...
        +  (G2_data.y_col(j) - G1_data.y_q_comb_outer).^2 );

    S21(j, :) = 1i.*G1_data.w_comb_outer.*besselh(0, k*dist(:))./4;

    % S22 computation
    t_lower = G2_data.t_grid_comb_outer(1:end - 1);
    t_upper = G2_data.t_grid_comb_outer(2:end);

    S22(j, :) = m1(k, G2_data.col_points(j), G2_data.t_mid_q_comb_outer, C1, C2)...
        .*w1_weights(k, G2_data.col_points(j), t_lower, t_upper) ...
        + G2_data.w_comb_outer.*m2(k, G2_data.col_points(j), ...
        G2_data.t_mid_q_comb_outer, C1, C2 );



end


% Now for the inner matricies that are more computationally expensive
K21 = zeros(length(G2_data.x_q_comb_outer), length(G1_data.x_q_comb_inner));

for j = 1:length(G2_data.x_q_comb_outer)

    dist = sqrt( (G2_data.x_q_comb_outer(j) - G1_data.x_q_comb_inner).^2 ...
        + (G2_data.y_q_comb_outer(j) - G1_data.y_q_comb_inner).^2);

    dist_dot_n = (G2_data.x_q_comb_outer(j) - G1_data.x_q_comb_inner)*G2_data.n(1) ...
        + (G2_data.y_q_comb_outer(j) - G1_data.y_q_comb_inner)*G2_data.n(2);

    K21(j, :) = (1i*k.*G1_data.w_comb_inner.*besselh(1, 1, k*dist).*dist_dot_n./dist)/2;


end

K21_inner = zeros(length(G2_data.x_q_comb_inner), length(G1_data.x_q_comb_inner));

for j = 1:length(G2_data.x_q_comb_inner)

    dist = sqrt( (G2_data.x_q_comb_inner(j) - G1_data.x_q_comb_inner).^2 ...
        + (G2_data.y_q_comb_inner(j) - G1_data.y_q_comb_inner).^2);

    dist_dot_n = (G2_data.x_q_comb_inner(j) - G1_data.x_q_comb_inner)*G2_data.n(1) ...
        + (G2_data.y_q_comb_inner(j) - G1_data.y_q_comb_inner)*G2_data.n(2);

    K21_inner(j, :) = (1i*k.*G1_data.w_comb_inner.*besselh(1, 1, k*dist).*dist_dot_n./dist)/2;


end

% Now for the inner matricies that are more computationally expensive
K12 = zeros(length(G1_data.x_q_comb_outer), length(G2_data.x_q_comb_inner));

for j = 1:length(G1_data.x_q_comb_outer)

    dist = sqrt( (G1_data.x_q_comb_outer(j) - G2_data.x_q_comb_inner).^2 ...
        + (G1_data.y_q_comb_outer(j) - G2_data.y_q_comb_inner).^2);

    dist_dot_n = (G1_data.x_q_comb_outer(j) - G2_data.x_q_comb_inner)*G1_data.n(1) ...
        + (G1_data.y_q_comb_outer(j) - G2_data.y_q_comb_inner)*G1_data.n(2);

    K12(j, :) = (1i*k.*G2_data.w_comb_inner.*besselh(1, 1, k*dist).*dist_dot_n./dist)/2;


end

K12_inner = zeros(length(G1_data.x_q_comb_inner), length(G2_data.x_q_comb_inner));

for j = 1:length(G1_data.x_q_comb_inner)

    dist = sqrt( (G1_data.x_q_comb_inner(j) - G2_data.x_q_comb_inner).^2 ...
        + (G1_data.y_q_comb_inner(j) - G2_data.y_q_comb_inner).^2);

    dist_dot_n = (G1_data.x_q_comb_inner(j) - G2_data.x_q_comb_inner)*G1_data.n(1) ...
        + (G1_data.y_q_comb_inner(j) - G2_data.y_q_comb_inner)*G1_data.n(2);

    K12_inner(j, :) = (1i*k.*G2_data.w_comb_inner.*besselh(1, 1, k*dist).*dist_dot_n./dist)/2;


end