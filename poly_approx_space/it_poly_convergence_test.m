% convergence testing of iterative method with comparison to all in 1
% polynomial solve

clear all
addpath('../General_functions/')

N_bf = [10, 20, 40, 80, 160, 320, 640] %, 1280];
R_max = 10;

% loading in true solution
all_in_one_solver = open('poly_solver_k10_Nbf_10_1280.mat')

% introducing screens
G1_data.G = all_in_one_solver.G1_data.G;
% G1_data.G = [-2*pi, 2*pi, 0, 0];

G2_data.G = all_in_one_solver.G2_data.G;
% G2_data.G = [2*pi, 0, 5*pi, 3*pi]; 

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


% % loading in true solution
% all_in_one_solver = open('poly_solver_k10_Nbf_10_1280.mat')
% 
% % introducing screens
% G1_data.G = all_in_one_solver.G1_data.G;
% % G1_data.G = [-2*pi, 2*pi, 0, 0];
% 
% G2_data.G = all_in_one_solver.G2_data.G;
% % G2_data.G = [2*pi, 0, 5*pi, 3*pi]; 
% 
% Lgrad_coeff = 0.15;
% alpha = 4;
% 
% % C_wl= 1/20; % degrees of freedom per wavelength for integration
% 
% k = 10;  % wavenumber
% 
% theta = 0;
% 
% % constants needed for the smoothing function
% C1 = 1;
% C2 = pi;


% pull out true solution we are measuring against

aj_1_true = all_in_one_solver.aj_1{end - 1};
aj_2_true = all_in_one_solver.aj_2{end - 1};

C_wl_err = 1/40;

[x1_1_q_err, y1_1_q_err, x1_2_q_err, y1_2_q_err, t1_grid_err, ...
    t1_mid_err, w1_err, N1_err, L1] = ...
    discretistion_vars_graded(G1_data.G, C_wl_err, k, Lgrad_coeff, alpha);

[x2_1_q_err, y2_1_q_err, x2_2_q_err, y2_2_q_err, t2_grid_err, ...
    t2_mid_err, w2_err, N2_err, L2] = ...
    discretistion_vars_graded(G2_data.G, C_wl_err, k, Lgrad_coeff, alpha);

t1_mid_err = t1_mid_err(100:end);
t2_mid_err = t2_mid_err(100:end);

x1_plotting = [t1_mid_err; flip(L1 - t1_mid_err) ];
x2_plotting = [t2_mid_err; flip(L2 - t2_mid_err) ];

w1_err = [w1_err(100:end); flip(w1_err(100:end))];
w2_err = [w2_err(100:end); flip(w2_err(100:end))];

% phi_1 = zeros(length(N_bf), length(t1_mid_err)*2);
% phi_2 = zeros(length(N_bf), length(t2_mid_err)*2);

phi_1_true = graded_coeff_2_solution(aj_1_true, ...
    all_in_one_solver.t1_bf_grid_store{end}, t1_mid_err, G1_data.L);

phi_2_true = graded_coeff_2_solution(aj_2_true, ...
    all_in_one_solver.t2_bf_grid_store{end}, t2_mid_err, G2_data.L);

phi_1_true_only_norm = sum(w1_err.*abs(phi_1_true));
phi_2_true_only_norm = sum(w2_err.*abs(phi_2_true));

% only adding one path to centralise solvers
% addpath('../General_functions/')

% % introducing screens
% G1_data.G = all_in_one_solver.G1_data.G;
% % G1_data.G = [-2*pi, 2*pi, 0, 0];
% 
% G2_data.G = all_in_one_solver.G2_data.G;
% % G2_data.G = [2*pi, 0, 5*pi, 3*pi]; 
% 
% Lgrad_coeff = 0.15;
% alpha = 4;

% C_wl= 1/20; % degrees of freedom per wavelength for integration
% 
% % k = 10;  % wavenumber
% % 
% % theta = 0;
% 
% % constants needed for the smoothing function
% C1 = 1;
% C2 = pi;
% 
% % quadrature nodes and other information needed
% [G1_data.x_1_q, G1_data.y_1_q, G1_data.x_2_q, G1_data.y_2_q, ...
%     G1_data.t_grid, G1_data.t_mid_q, G1_data.w, G1_data.N, G1_data.L] = ...
%     discretistion_vars_graded(G1_data.G, C_wl, k, Lgrad_coeff, alpha);
% 
% [G2_data.x_1_q, G2_data.y_1_q, G2_data.x_2_q, G2_data.y_2_q, ...
%     G2_data.t_grid, G2_data.t_mid_q, G2_data.w, G2_data.N, G2_data.L] = ...
%     discretistion_vars_graded(G2_data.G, C_wl, k, Lgrad_coeff, alpha);


% N_bf = [10, 20, 40, 80, 160] %, 320, 640, 1280];
% R_max = 20;

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
    
    
    [aj_1_R{j}, aj_2_R{j}, phi_1_r{j}, phi_2_r{j}] = iterative_poly_graded_PIM_solve(...
        S11, S12, S21, S22, u_inc1, u_inc2, R_max, t1_bf_grid, t2_bf_grid,...
        G1_data, G2_data);
    toc

    % Error computation

    % First computing approximation at the approximation for each r and set
    % of basis fucntions
    for r = 1:R_max

        phi_1_r{j, r} = graded_coeff_2_solution(aj_1_R{j}(:, r), t1_bf_grid, ...
            t1_mid_err, G1_data.L);
        
        phi_2_r{j, r} = graded_coeff_2_solution(aj_2_R{j}(:, r), t2_bf_grid, ...
            t2_mid_err, G2_data.L);

        err_1_r{j}(r) = sum(w1_err.*abs(phi_1_r{j, r} - phi_1_true))./phi_1_true_only_norm;

        err_2_r{j}(r) = sum(w2_err.*abs(phi_2_r{j, r} - phi_2_true))./phi_2_true_only_norm;

    end

    t1_bf_grid_store{j} = t1_bf_grid;
    t2_bf_grid_store{j} = t2_bf_grid;


end


%%
% plotting the error, we want to show two things, error decreasing as a
% function of N_bf and also as a function of R
figure()
r1 = [0:2:R_max*2-2];
r2 = [1:2:2*R_max - 1];
for n = 1:length(N_bf)
    txt = 2*length(t1_bf_grid_store{n});
    txt1 = "N = " + txt;
    subplot(2, 1, 1)
    plot(r1, err_1_r{n}, 'DisplayName', txt1)
    hold on
    xlabel('r')
    ylabel('Error')
    legend show

    txt2 = "N = " + 2*length(t2_bf_grid_store{n});
    subplot(2, 1, 2)
    plot(r2, err_2_r{n}, 'DisplayName', txt2)
    hold on
    xlabel('r')
    ylabel('Error')
    legend show

end


