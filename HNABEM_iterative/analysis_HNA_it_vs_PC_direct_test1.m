% comparison between HNA iterative and piecewise constant standard BEM
% solver

clear all

% adding paths

addpath('../General_functions/')
addpath('../../BEAM_HNABEMLAB/')
addPathsHNA  % allows HNABEM to find all of the relevatn subfolders

% loading in data from HNA solver

load('test3_HNA_pmax4_overlap2.mat')

G1_data.G = G1_data_HNA{1}.G;
G2_data.G = G2_data_HNA{1}.G;


%% loading in data from piecewise constant direct solver
 
addpath('../poly_approx_space/')
load('PC_direct_test3_k10_theta0.mat')

G1_data_PCD.G = info_needed.G1;
G2_data_PCD.G = info_needed.G2;

k = info_needed.k;
theta = info_needed.theta;
Lgrad_coeff_poly = info_needed.L_grad_coeff;
alpha_poly = info_needed.alpha;
bf_dof_per_wl = info_needed.bf_dof_per_wl;

G1_data_PCD = get_bf_graded_grid(G1_data_PCD, bf_dof_per_wl(end), k, ...
    Lgrad_coeff_poly, alpha_poly);

G2_data_PCD = get_bf_graded_grid(G2_data_PCD, bf_dof_per_wl(end), k, ...
    Lgrad_coeff_poly, alpha_poly);

d = [sin(theta) -cos(theta) ];

%%%%% For this to be a good comparison we need to choose at least as many
%%%%% points as either the HNA or PC basis functions. Now we know that the
%%%%% PC is going to be larger. Seen as we are generally using 80 dof per
%%%%% wavelength for the true quadrature solution, we should probably use
%%%%% 100-120 dof per wavelength for the mesh for the error computations.
%%%%% Then we need to compute PC phi and HNA phi at those points for the
%%%%% maximum iteration.
%% Creating grid
eval_at_and_err_dof_per_wl = 1/20;
int_for_sol_dof_per_wl_outer = 1/20;
int_for_sol_dof_per_wl_inner = 1/20;


G1_data_err.G = G1_data.G;
G1_data_err.L = sqrt( (G1_data_err.G(3) - G1_data_err.G(1))^2 +...
    (G1_data_err.G(4) - G1_data_err.G(2))^2 );

G2_data_err.G = G2_data.G;
G2_data_err.L = sqrt( (G2_data_err.G(3) - G2_data_err.G(1))^2 +...
    (G2_data_err.G(4) - G2_data_err.G(2))^2 );

% G1_data_err = get_bf_graded_grid(G1_data_err, err_dof_per_wl, k, ...
%     Lgrad_coeff_poly, alpha_poly);
% 
% G2_data_err = get_bf_graded_grid(G2_data_err, err_dof_per_wl, k, ...
%     Lgrad_coeff_poly, alpha_poly);

% extra bits for HNA part needed
G1_data_err.n = [-(G1_data_err.G(4) - G1_data_err.G(2)), G1_data_err.G(3) - G1_data_err.G(1)]...
    /G1_data_err.L;
G2_data_err.n = [-(G2_data_err.G(4) - G2_data_err.G(2)), G2_data_err.G(3) - G2_data_err.G(1)]...
    /G1_data_err.L;

G1_data_err.alpha = -sign(dot(d, G1_data_err.n));
G2_data_err.alpha = -sign(dot(d, G2_data_err.n));

G1_data_int = G1_data_err;
G2_data_int = G2_data_err;

% Quadrature - for computing solution at certain points and for errors
% integration
G1_data_err = get_graded_quad_points_HF_it(G1_data_err, eval_at_and_err_dof_per_wl,...
    1/10, k, Lgrad_coeff_poly, alpha_poly);

G2_data_err = get_graded_quad_points_HF_it(G2_data_err, eval_at_and_err_dof_per_wl,...
    1/10, k, Lgrad_coeff_poly, alpha_poly);

% Quadrature - for computing phi1 and phi2 at the points computed above
G1_data_int = get_graded_quad_points_HF_it(G1_data_int, int_for_sol_dof_per_wl_outer,...
    int_for_sol_dof_per_wl_inner, k, Lgrad_coeff_poly, alpha_poly);

G2_data_int = get_graded_quad_points_HF_it(G2_data_int, int_for_sol_dof_per_wl_outer,...
    int_for_sol_dof_per_wl_inner, k, Lgrad_coeff_poly, alpha_poly);


%% compute whether on the positive or negative side of the screen

G2_data_int = get_coeff_for_different_sides_of_screen(G1_data_int, ...
    G2_data_int);
G1_data_int = get_coeff_for_different_sides_of_screen(G2_data_int, ...
    G1_data_int);


%% compute HNA solution


    
%% functionising
for n = 1:length(bf_dof_per_wl)
    tic
    [phi1, phi2] = get_HNA_it_soln_at_points(G1_data_int, G2_data_int, k, ...
        d, v_N1_HNA_cell{n}, v_N2_HNA_cell{n}, length(phi1_HNA{1}));
    toc
    phi1_HNA_eval(:, n) = phi1(:, end);
    phi2_HNA_eval(:, n) = phi2(:, end);
end

%% compute PC solution


phi1_PCD =  graded_coeff_2_solution(aj1_coeff{end}, ...
    G1_data_PCD.t_bf_grid, G1_data_int.t_mid_q_outer, ...
    G1_data_PCD.L);

phi2_PCD =  graded_coeff_2_solution(aj2_coeff{end}, ...
    G2_data_PCD.t_bf_grid, G2_data_int.t_mid_q_outer, ...
    G2_data_PCD.L);

%% compute the error

for n = 1:length(bf_dof_per_wl)
    err_phi1(n) = sum((abs(phi1_PCD - phi1_HNA_eval(:, n))./abs(phi1_PCD))...
        .*G1_data_int.w_comb_outer);

    err_phi1_test(n) = sum(abs(phi1_PCD - phi1_HNA_eval(:, n)).*G1_data_int.w_comb_outer)./sum(abs(phi1_PCD)...
        .*G1_data_int.w_comb_outer);
    
    err_phi2(n) = sum((abs(phi2_PCD - phi2_HNA_eval(:, n))./abs(phi2_PCD))...
        .*G2_data_int.w_comb_outer);

end

err_phi1

err_phi2
%%
figure()
plot(G1_data_int.t_mid_q_comb_outer/G1_data_int.L, real(phi1_PCD), 'DisplayName', 'Direct')
hold on
plot(G1_data_int.t_mid_q_comb_outer/G1_data_int.L, real(phi1_HNA_eval(:, end)), 'DisplayName', 'HNA')
legend show

figure()
plot(G2_data_int.t_mid_q_comb_outer/G2_data_int.L, real(phi2_PCD), 'DisplayName', 'Direct')
hold on
plot(G2_data_int.t_mid_q_comb_outer/G2_data_int.L, real(phi2_HNA_eval(:, end)), 'DisplayName', 'HNA')
legend show

%% Old:
% phi1_HNA_eval = zeros(length(G1_data_err.t_mid_q_comb_outer), length(phi1_HNA{1}));
% phi2_HNA_eval = zeros(length(G2_data_err.t_mid_q_comb_outer), length(phi2_HNA{1}));
% 
% phi1_HNA_eval(:, 1) = v_N1_HNA_cell{1}{1}.eval(G1_data_err.t_mid_q_comb_outer, 1) ...
%     + 2*G1_data_err.alpha*duidn(G1_data_err, G1_data_err.L, k, d, ...
%     G1_data_err.t_mid_q_comb_outer);
% 
% phi1_HNA_eval_inner = zeros(length(G1_data_int.t_mid_q_comb_inner), length(phi1_HNA{1}));
% phi2_HNA_eval_inner = zeros(length(G2_data_int.t_mid_q_comb_inner), length(phi2_HNA{1}));
% 
% phi1_HNA_eval_inner(:, 1) = v_N1_HNA_cell{1}{1}.eval(G1_data_int.t_mid_q_comb_inner, 1) ...
%     + 2*G1_data_int.alpha*duidn(G1_data_int, G1_data_int.L, k, d, ...
%     G1_data_int.t_mid_q_comb_inner);
% 
% for r = 2:length(phi1_HNA{1})
% 
%     % compute phi2 outer
%     phi2_HNA_eval(:, r - 1) = v_N2_HNA_cell{1}{r-1}.eval(G2_data_err.t_mid_q_comb_outer, 1) +...
%         2*G2_data_err.alpha*duidn(G2_data_err, G2_data_err.L, k, d, ...
%         G2_data_err.t_mid_q_comb_outer) ...
%         + midpoint_dphikdn_f_diff_screen(k, ...
%         G2_data_err.x_q_comb_outer, G2_data_err.y_q_comb_outer, ...
%         G1_data_int.w_comb_inner, G1_data_int.x_q_comb_inner, ...
%         G1_data_int.y_q_comb_inner, G1_data_int.beta_inner.*phi1_HNA_eval_inner(:, r - 1), ...
%         G2_data_err.n);
% 
%     % compute phi2 inner
%     phi2_HNA_eval_inner(:, r - 1) = v_N2_HNA_cell{1}{r-1}.eval(G2_data_int.t_mid_q_comb_inner, 1) +...
%         2*G2_data_int.alpha*duidn(G2_data_int, G2_data_int.L, k, d, ...
%         G2_data_int.t_mid_q_comb_inner) ...
%         + midpoint_dphikdn_f_diff_screen(k, ...
%         G2_data_int.x_q_comb_inner, G2_data_int.y_q_comb_inner, ...
%         G1_data_int.w_comb_inner, G1_data_int.x_q_comb_inner, ...
%         G1_data_int.y_q_comb_inner, G1_data_int.beta_inner.*phi1_HNA_eval_inner(:, r-1), ...
%         G2_data_int.n);
% 
%     % computing the phi1 outer
%  phi1_HNA_eval(:, r) = v_N1_HNA_cell{1}{r}.eval(G1_data_err.t_mid_q_comb_outer, 1) +...
%         2*G1_data_err.alpha*duidn(G1_data_err, G1_data_err.L, k, d, ...
%         G1_data_err.t_mid_q_comb_outer) ...
%         + midpoint_dphikdn_f_diff_screen(k, ...
%         G1_data_err.x_q_comb_outer, G1_data_err.y_q_comb_outer, ...
%         G2_data_int.w_comb_inner, G2_data_int.x_q_comb_inner, ...
%         G2_data_int.y_q_comb_inner, G2_data_int.beta_inner.*phi2_HNA_eval_inner(:, r - 1), ...
%         G1_data_err.n);
% 
%  % phi1 inner
%      phi1_HNA_eval_inner(:, r) = v_N1_HNA_cell{1}{r}.eval(G1_data_int.t_mid_q_comb_inner, 1) +...
%         2*G1_data_int.alpha*duidn(G1_data_int, G1_data_int.L, k, d, ...
%         G1_data_int.t_mid_q_comb_inner) ...
%         + midpoint_dphikdn_f_diff_screen(k, ...
%         G1_data_int.x_q_comb_inner, G1_data_int.y_q_comb_inner, ...
%         G2_data_int.w_comb_inner, G2_data_int.x_q_comb_inner, ...
%         G2_data_int.y_q_comb_inner, G2_data_int.beta_inner.*phi2_HNA_eval_inner(:, r - 1), ...
%         G1_data_int.n);
% 
% 
% 
%     
% 
% end
% 
% % final phi2 calc
% 
% phi2_HNA_eval(:, end) = v_N2_HNA_cell{1}{end}.eval(G2_data_err.t_mid_q_comb_outer, 1) +...
%         2*G2_data_err.alpha*duidn(G2_data_err, G2_data_err.L, k, d, ...
%         G2_data_err.t_mid_q_comb_outer) ...
%         + midpoint_dphikdn_f_diff_screen(k, ...
%         G2_data_err.x_q_comb_outer, G2_data_err.y_q_comb_outer, ...
%         G1_data_int.w_comb_inner, G1_data_int.x_q_comb_inner, ...
%         G1_data_int.y_q_comb_inner, G1_data_int.beta_inner.*phi1_HNA_eval_inner(:, end), ...
%         G2_data_err.n);




