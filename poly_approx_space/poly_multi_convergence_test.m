% polynomial all in one convergence test

clear all
tic
% only adding one path to centralise solvers
addpath('../General_functions/')

% switch and things get interesting
G1_data.G = [-2*pi, 2*pi, 0, 0];
% G1_data.L = sqrt( (G1(3) - G1(1))^2 +(G1(4) - G1(2))^2 );

G2_data.G = [2*pi, 0, 5*pi, 3*pi]; 
% G2_data.L = sqrt( (G2(3) - G2(1))^2 +(G2(4) - G2(2))^2 );

Lgrad_coeff = 0.15;
alpha = 4;

C_wl= 1/20;

k = 10;  % wavenumber

theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

% run_max = 5;

% quadrature nodes and other information needed
[G1_data.x_1_q, G1_data.y_1_q, G1_data.x_2_q, G1_data.y_2_q, G1_data.t_grid, G1_data.t_mid, G1_data.w, G1_data.N, G1_data.L] = ...
    discretistion_vars_graded(G1_data.G, C_wl, k, Lgrad_coeff, alpha);

[G2_data.x_1_q, G2_data.y_1_q, G2_data.x_2_q, G2_data.y_2_q, G2_data.t_grid, G2_data.t_mid, G2_data.w, G2_data.N, G2_data.L] = ...
    discretistion_vars_graded(G2_data.G, C_wl, k, Lgrad_coeff, alpha);

% G1_data.s = [ G1_data.t_mid(1:end) ; flip(G1_data.L - G1_data.t_mid(1:end)) ];
% G1_data.x_col = [ G1_data.x_1_q(1:end) ; flip(G1_data.x_2_q(1:end)) ];
% G1_data.y_col = [ G1_data.y_1_q(1:end) ; flip(G1_data.y_2_q(1:end)) ];
% 
% % col_choice2 = sort(randi(length(G2_data.t_mid(:)), 20, 1));
% G2_data.s = [ G2_data.t_mid(1:end) ; flip(G2_data.L - G2_data.t_mid(1:end)) ];
% G2_data.x_col = [ G2_data.x_1_q(1:end) ; flip(G2_data.x_2_q(1:end)) ];
% G2_data.y_col = [ G2_data.y_1_q(1:end) ; flip(G2_data.y_2_q(1:end)) ];

N_bf = [10, 20, 40, 80, 160, 320, 640] %, 1280, 2560];

for j = 1:length(N_bf)
    
    disp(j)
    % The support for the basis functions
%     C_wl_bf1 = 1/(j+1);
%     C_wl_bf2 = 1/(j+1);
    [x1_1, y1_1, x1_2, y1_2, t1_bf_grid, t1_mid, ~, ~] = ...
    given_N_discretistion_vars_graded(G1_data.G, N_bf(j), Lgrad_coeff, alpha);
    [x2_1, y2_1, x2_2, y2_2, t2_bf_grid, t2_mid, ~, ~] = ...
    given_N_discretistion_vars_graded(G2_data.G, N_bf(j), Lgrad_coeff, alpha);

    % collocation points
    G1_data.s = [ t1_mid ; flip(G1_data.L - t1_mid) ];
    G1_data.x_col = [ x1_1 ; flip(x1_2) ];
    G1_data.y_col = [ y1_1 ; flip(y1_2) ];
    
    % col_choice2 = sort(randi(length(G2_data.t_mid(:)), 20, 1));
    G2_data.s = [ t2_mid ; flip(G2_data.L - t2_mid) ];
    G2_data.x_col = [ x2_1 ; flip(x2_2) ];
    G2_data.y_col = [ y2_1 ; flip(y2_2) ];
        
    while length(G1_data.t_grid) < length(t1_bf_grid)  % case where there are more basis function than integration points, rendering the new bf useless
        C_wl= C_wl/2;
        % quadrature nodes and other information needed
    [G1_data.x_1_q, G1_data.y_1_q, G1_data.x_2_q, G1_data.y_2_q, G1_data.t_grid, G1_data.t_mid, G1_data.w, G1_data.N, G1_data.L] = ...
        discretistion_vars_graded(G1_data.G, C_wl, k, Lgrad_coeff, alpha);
    
    [G2_data.x_1_q, G2_data.y_1_q, G2_data.x_2_q, G2_data.y_2_q, G2_data.t_grid, G2_data.t_mid, G2_data.w, G2_data.N, G2_data.L] = ...
        discretistion_vars_graded(G2_data.G, C_wl, k, Lgrad_coeff, alpha);
    
%     G1_data.s = [ G1_data.t_mid(1:end) ; flip(G1_data.L - G1_data.t_mid(1:end)) ];
%     G1_data.x_col = [ G1_data.x_1_q(1:end) ; flip(G1_data.x_2_q(1:end)) ];
%     G1_data.y_col = [ G1_data.y_1_q(1:end) ; flip(G1_data.y_2_q(1:end)) ];
%     
%     
%     G2_data.s = [ G2_data.t_mid(1:end) ; flip(G2_data.L - G2_data.t_mid(1:end)) ];
%     G2_data.x_col = [ G2_data.x_1_q(1:end) ; flip(G2_data.x_2_q(1:end)) ];
%     G2_data.y_col = [ G2_data.y_1_q(1:end) ; flip(G2_data.y_2_q(1:end)) ];
    disp('Increasing number of quadrature points due to basis functions')
    end
    [S11, S12, S21, S22, u_inc1, u_inc2] = ...
        compute_matrices_for_iterative_solve(G1_data, G2_data, k, ...
        t1_bf_grid, t2_bf_grid, theta, C1, C2 );
    
    % all in one solve
    A = [S11  S12 ; S21 S22];
    
    u_inc = [u_inc1 ; u_inc2];
    
    coeffs = A\u_inc;

    t1_bf_grid_store{j, :} = t1_bf_grid;
    t2_bf_grid_store{j, :} = t2_bf_grid;
    
    aj_1{j, :} = coeffs(1:2*length(t1_bf_grid) - 2);
    aj_2{j, :} = coeffs(2*length(t1_bf_grid) - 2+ 1: end);

    clear A S11 S12 S21 S22 uinc u_inc1 u_inc2

end


%%
% x1_plot = linspace(0.01, G1_data.L/2 - 0.01, 300);
% w1_plot = (x1_plot(end) - x1_plot(1))/(300+1);
% 
% x1_plotting = [x1_plot(:) ; (G1_data.L/2 + x1_plot(:) )];
% 
% x2_plot = linspace(0.01, G2_data.L/2 - 0.01, 300);
% w2_plot = (x2_plot(end) - x2_plot(1))/(300+1);
% 
% x2_plotting = [x2_plot(:) ; (G2_data.L/2 + x2_plot(:) )];

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

phi_1 = zeros(length(N_bf), length(t1_mid_err)*2);
phi_2 = zeros(length(N_bf), length(t2_mid_err)*2);

for j = 1:length(N_bf)

    phi_1(j, :) = graded_coeff_2_solution(aj_1{j, :}, t1_bf_grid_store{j},...
        t1_mid_err, G1_data.L);
    phi_2(j, :) =graded_coeff_2_solution(aj_2{j, :}, t2_bf_grid_store{j},...
        t2_mid_err, G2_data.L);

end

% plotting

figure()
for j = 1:length(N_bf)
    plot(x1_plotting/G1_data.L, phi_1(j,:), 'DisplayName', ...
        strcat('N = ', num2str(length(aj_1{j}))))
    hold on

end
xlabel('$x/L_1$')
ylabel('$\phi_1(x)$')
title('Polynoial approximation solve, $\phi_{1}(x)$')
legend show
xlim([-0.05 1.05])

figure()
for j = 1:length(N_bf)
    plot(x2_plotting/G2_data.L, phi_2(j,:), 'DisplayName', ...
        strcat('N = ', num2str(length(aj_2{j}))))
    hold on

end
xlabel('$x/L_2$')
ylabel('$\phi_{2}(x)$')
title('Polynoial approximation solve, $\phi_{2}(x)$')
legend show
xlim([-0.05 1.05])

%% Error/ convergence plot
phi_1_only_norm = sum(w1_err.'.*abs(phi_1(end, :)));
phi_2_only_norm = sum(w2_err.'.*abs(phi_2(end, :)));

for j = 1:length(N_bf) - 1
    
    err_phi_1(j) = sum(w1_err.'.*abs(phi_1(end, :) - phi_1(j, :)))./phi_1_only_norm;
        

    err_phi_2(j) = sum(w2_err.'.*abs(phi_2(end, :) - phi_2(j, :)))./phi_2_only_norm;

end

%% EOC

for j = 1:length(N_bf) - 2

    EOC_phi_1(j) = log2(err_phi_1(j)/err_phi_1(j+1));

    EOC_phi_2(j) = log2(err_phi_2(j)/err_phi_2(j+1));

end

err_phi_1, err_phi_2, EOC_phi_1, EOC_phi_2

toc
% 
filename = 'poly_solver_k10_Nbf_10_640.mat';

save(filename, "G1_data", "G2_data", "N_bf", "aj_1", "aj_2", "t1_bf_grid_store", "t2_bf_grid_store")
% % 
