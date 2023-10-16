% convergence testing of iterative method with comparison to all in 1
% polynomial solve

clear all

% only adding one path to centralise solvers
addpath('../General_functions/')

% introducing screens
G1_data.G = [-2*pi, 2*pi, 0, 0];

G2_data.G = [2*pi, 0, 5*pi, 3*pi]; 

Lgrad_coeff = 0.15;
alpha = 4;

C_wl= 1/20; % degrees of freedom per wavelength for integration

k = 10;  % wavenumber

theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

% quadrature nodes and other information needed
[G1_data.x_1_q, G1_data.y_1_q, G1_data.x_2_q, G1_data.y_2_q, ...
    G1_data.t_grid, G1_data.t_mid_q, G1_data.w, G1_data.N, G1_data.L] = ...
    discretistion_vars_graded(G1_data.G, C_wl, k, Lgrad_coeff, alpha);

[G2_data.x_1_q, G2_data.y_1_q, G2_data.x_2_q, G2_data.y_2_q, ...
    G2_data.t_grid, G2_data.t_mid_q, G2_data.w, G2_data.N, G2_data.L] = ...
    discretistion_vars_graded(G2_data.G, C_wl, k, Lgrad_coeff, alpha);


N_bf = [10, 20, 40, 80, 160, 320, 640, 1280];

for j = 1:length(N_bf)
    tic
    disp(j)
    [G1_data.x_1_col, G1_data.y_1_col, G1_data.x_2_col, G1_data.y_2_col,...
        t1_bf_grid, G1_data.t_mid_col, ~, ~] =  ...
        given_N_discretistion_vars_graded(G1_data.G, N_bf(j), Lgrad_coeff, alpha);
    
    [G2_data.x_1_col, G2_data.y_1_col, G2_data.x_2_col, G2_data.y_2_col,...
        t2_bf_grid, G2_data.t_mid_col, ~, ~] = ...
        given_N_discretistion_vars_graded(G2_data.G, N_bf(j), Lgrad_coeff, alpha);
  
    G1_data.s = [ G1_data.t_mid_col(1:end) ; flip(G1_data.L - ...
        G1_data.t_mid_col(1:end)) ];
    G1_data.x_col = [ G1_data.x_1_col(1:end) ; flip(G1_data.x_2_col(1:end)) ];
    G1_data.y_col = [ G1_data.y_1_col(1:end) ; flip(G1_data.y_2_col(1:end)) ];
    
    % col_choice2 = sort(randi(length(G2_data.t_mid(:)), 20, 1));
    G2_data.s = [ G2_data.t_mid_col(1:end) ; flip(G2_data.L - ...
        G2_data.t_mid_col(1:end)) ];
    G2_data.x_col = [ G2_data.x_1_col(1:end) ; flip(G2_data.x_2_col(1:end)) ];
    G2_data.y_col = [ G2_data.y_1_col(1:end) ; flip(G2_data.y_2_col(1:end)) ];


%     [x1_1, y1_1, x1_2, y1_2, t1_bf_grid, t1_mid, ~, ~] = ...
%     given_N_discretistion_vars_graded(G1_data.G, N_bf(j), Lgrad_coeff, alpha);
%     [x2_1, y2_1, x2_2, y2_2, t2_bf_grid, t2_mid, ~, ~] = ...
%     given_N_discretistion_vars_graded(G2_data.G, N_bf(j), Lgrad_coeff, alpha);
% 
%     % collocation points
%     G1_data.s = [ t1_mid ; flip(G1_data.L - t1_mid) ];
%     G1_data.x_col = [ x1_1 ; flip(x1_2) ];
%     G1_data.y_col = [ y1_1 ; flip(y1_2) ];
%     
%     % col_choice2 = sort(randi(length(G2_data.t_mid(:)), 20, 1));
%     G2_data.s = [ t2_mid ; flip(G2_data.L - t2_mid) ];
%     G2_data.x_col = [ x2_1 ; flip(x2_2) ];
%     G2_data.y_col = [ y2_1 ; flip(y2_2) ];
        
    while length(G1_data.t_grid) < length(t1_bf_grid)  % case where there are more basis function than integration points, rendering the new bf useless
        C_wl= C_wl/2;
        % quadrature nodes and other information needed
        [G1_data.x_1_q, G1_data.y_1_q, G1_data.x_2_q, G1_data.y_2_q, ...
            G1_data.t_grid, G1_data.t_mid_q, G1_data.w, G1_data.N, G1_data.L] = ...
            discretistion_vars_graded(G1_data.G, C_wl, k, Lgrad_coeff, alpha);
        
        [G2_data.x_1_q, G2_data.y_1_q, G2_data.x_2_q, G2_data.y_2_q, ...
            G2_data.t_grid, G2_data.t_mid_q, G2_data.w, G2_data.N, G2_data.L] = ...
            discretistion_vars_graded(G2_data.G, C_wl, k, Lgrad_coeff, alpha);
          
        disp('Increasing number of quadrature points due to basis functions')
    end
    [S11, S12, S21, S22, u_inc1, u_inc2] = ...
        compute_matrices_for_iterative_solve(G1_data, G2_data, k, ...
        t1_bf_grid, t2_bf_grid, theta, C1, C2 );

    % iterative solve
    R_max = 20;
    
    [aj_1_R{j}, aj_2_R{j}, phi_1_r{j}, phi_2_r{j}] = iterative_poly_graded_PIM_solve(...
        S11, S12, S21, S22, u_inc1, u_inc2, R_max, t1_bf_grid, t2_bf_grid,...
        G1_data, G2_data);
    toc

end




