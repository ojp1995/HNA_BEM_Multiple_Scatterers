% Comparison between piecewise constant solvers 
%%%% TEST 5a
clear all

addpath('../General_functions/')

% Now we want to load each dataset, starting with the direct polynomial
% solver

load('PC_direct_test5b_k10_thetapi_4_short.mat')

% PCD - Piecewise constant direct solver
aj1_coeff_PCD = aj1_coeff;
aj2_coeff_PCD = aj2_coeff;

G1_data.G = info_needed.G1;
G2_data.G = info_needed.G2;

L_grad_coeff = info_needed.L_grad_coeff;
alpha = info_needed.alpha;

k = info_needed.k;

bf_dof_per_wl = info_needed.bf_dof_per_wl;

% Now load in iterative method
load('it_PC_direct_test5b_k10_thetapi_4_short.mat')

aj1_it_coeff = aj1_coeff;
aj2_it_coeff = aj2_coeff;

R_max = 150;

% Now we need to compute the points we want to evaluate phi at for each of
% them
G1_err_quad = get_bf_graded_grid(G1_data, bf_dof_per_wl(end)/2, k, ...
        L_grad_coeff, alpha);
G2_err_quad = get_bf_graded_grid(G2_data, bf_dof_per_wl(end)/2, k, ...
        L_grad_coeff, alpha);

[~, ~, ~, ~, ~, ~, G1_err_quad.w, ~, ~] = ...
    discretistion_vars_graded(G1_err_quad.G, bf_dof_per_wl(end)/2, k, ...
    L_grad_coeff, alpha);

[~, ~, ~, ~, ~, ~, G2_err_quad.w, ~, ~] = ...
    discretistion_vars_graded(G2_err_quad.G, bf_dof_per_wl(end)/2, k, ...
    L_grad_coeff, alpha);

% Now compute the grids for PCD solver (only most refined solution as that
% is what we are interested in)

G1_data_PCD = get_bf_graded_grid(G1_data, bf_dof_per_wl(end), k, ...
        L_grad_coeff, alpha);

G2_data_PCD = get_bf_graded_grid(G2_data, bf_dof_per_wl(end), k, ...
        L_grad_coeff, alpha);

x1_plotting = linspace(0.01, G1_data_PCD.L/2, 1000);
x2_plotting = linspace(0.01, G2_data_PCD.L/2, 1000);

% compute phi1 and phi2 at points of interest - for plotting
phi1_PCD_plotting = graded_coeff_2_solution(aj1_coeff_PCD{end}, ...
        G1_data_PCD.t_bf_grid, x1_plotting, ...
        G1_data_PCD.L);

phi2_PCD_plotting = graded_coeff_2_solution(aj2_coeff_PCD{end}, ...
        G2_data_PCD.t_bf_grid, x2_plotting, ...
        G2_data_PCD.L);

% compute phi1 and phi2 at points of interest - for error
phi1_PCD_for_err = graded_coeff_2_solution(aj1_coeff_PCD{end}, ...
        G1_data_PCD.t_bf_grid, G1_err_quad.t_mid_col, ...
        G1_data_PCD.L);

phi2_PCD_for_err = graded_coeff_2_solution(aj2_coeff_PCD{end}, ...
        G2_data_PCD.t_bf_grid, G2_err_quad.t_mid_col, ...
        G2_data_PCD.L);

%% sanity check
figure()
plot([x1_plotting (G1_data_PCD.L - flip(x1_plotting))]/G1_data_PCD.L, real(phi1_PCD_plotting),...
    'DisplayName', 'PC Direct 80 dof per wl')

figure()
plot([x2_plotting (G2_data_PCD.L - flip(x2_plotting))]/G2_data_PCD.L, real(phi2_PCD_plotting),...
    'DisplayName', 'PC Direct 80 dof per wl')

%% comparison to iterative solver
err_1 = zeros(length(bf_dof_per_wl), R_max);
err_2 = zeros(length(bf_dof_per_wl), R_max);

for n = 1:length(bf_dof_per_wl)  % looping over the basis functions first
    % create the meshes for this number of dof
    G1_data_it = get_bf_graded_grid(G1_data, bf_dof_per_wl(n), k, ...
        L_grad_coeff, alpha);

    G2_data_it = get_bf_graded_grid(G2_data, bf_dof_per_wl(n), k, ...
        L_grad_coeff, alpha);

%     phi1_PCD_for_err = graded_coeff_2_solution(aj1_coeff_PCD{n}, ...
%         G1_data_it.t_bf_grid, G1_err_quad.t_mid_col, ...
%         G1_data_PCD.L);
% 
%     phi2_PCD_for_err = graded_coeff_2_solution(aj2_coeff_PCD{n}, ...
%         G2_data_it.t_bf_grid, G2_err_quad.t_mid_col, ...
%         G2_data_PCD.L);


    for r = R_max:R_max  % now computing the value of phi at the quadrature points to then integrate

        phi_1_it = graded_coeff_2_solution(aj1_it_coeff{n}(:, r), ...
        G1_data_it.t_bf_grid, G1_err_quad.t_mid_col, ...
        G1_data_it.L);

        phi_2_it = graded_coeff_2_solution(aj2_it_coeff{n}(:, r), ...
        G2_data_it.t_bf_grid, G2_err_quad.t_mid_col, ...
        G2_data_it.L);

        % now compute the error:
        err_1(n, r) = sum( (abs(phi_1_it - phi1_PCD_for_err)...
            ./abs(phi1_PCD_for_err)).*[G1_err_quad.w ;...
            flip(G1_err_quad.w)] );

         
        err_2(n, r) = sum( (abs(phi_2_it - phi2_PCD_for_err)...
            ./abs(phi2_PCD_for_err)).*[G2_err_quad.w ;...
            flip(G2_err_quad.w)] );

    end

end

