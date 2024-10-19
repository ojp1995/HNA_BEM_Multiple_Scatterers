% comparison between direct method and higher order solver

clear all

load('PC_direct_test1_k10_theta0.mat')

k = info_needed.k;
theta = info_needed.theta;
Lgrad_coeff_poly = info_needed.L_grad_coeff;
alpha_poly = info_needed.alpha;
bf_dof_per_wl = info_needed.bf_dof_per_wl;

G1_data_poly.G = info_needed.G1;
G2_data_poly.G = info_needed.G2;

G1_data_poly_true =  get_bf_graded_grid(G1_data_poly, bf_dof_per_wl(end), k, ...
        Lgrad_coeff_poly, alpha_poly);
                
G2_data_poly_true =  get_bf_graded_grid(G2_data_poly, bf_dof_per_wl(end), k, ...
        Lgrad_coeff_poly, alpha_poly);

phi_1_poly_true = graded_coeff_2_solution(aj1_coeff{end}, ...
        G1_data_poly_true.t_bf_grid, G1_data_poly_true.t_mid_col, ...
        G1_data_poly_true.L);

phi_2_poly_true = graded_coeff_2_solution(aj2_coeff{end}, ...
        G2_data_poly_true.t_bf_grid, G2_data_poly_true.t_mid_col, ...
        G2_data_poly_true.L);

% 10 dof solution

G1_data_poly_10dof =  get_bf_graded_grid(G1_data_poly, bf_dof_per_wl(2), k, ...
        Lgrad_coeff_poly, alpha_poly);
                
G2_data_poly_10dof =  get_bf_graded_grid(G2_data_poly, bf_dof_per_wl(2), k, ...
        Lgrad_coeff_poly, alpha_poly);

phi_1_poly_10dof = graded_coeff_2_solution(aj1_coeff{2}, ...
        G1_data_poly_10dof.t_bf_grid, G1_data_poly_true.t_mid_col, ...
        G1_data_poly_10dof.L);

phi_2_poly_10dof = graded_coeff_2_solution(aj2_coeff{2}, ...
        G2_data_poly_10dof.t_bf_grid, G2_data_poly_true.t_mid_col, ...
        G2_data_poly_10dof.L);

%%
x1 = [G1_data_poly_true.t_mid_col ;flip(G1_data_poly_true.L - G1_data_poly_true.t_mid_col)];
x2 = [G2_data_poly_true.t_mid_col ;flip(G2_data_poly_true.L - G2_data_poly_true.t_mid_col)];

figure();
plot(x1, real(phi_1_poly_10dof), 'DisplayName', '10 dof')
hold on
plot(x1, real(phi_1_poly_true), 'DisplayName', '80 dof')
legend show

figure();
plot(x2, real(phi_2_poly_10dof), 'DisplayName', '10 dof')
hold on
plot(x2, real(phi_2_poly_true), 'DisplayName', '80 dof')
legend show


%% error computation
weights1 = [G1_data_poly_true.t_bf_grid(2:end) - G1_data_poly_true.t_bf_grid(1:end-1);...
    flip(G1_data_poly_true.t_bf_grid(2:end) - G1_data_poly_true.t_bf_grid(1:end-1))];

weights2 = [G2_data_poly_true.t_bf_grid(2:end) - G2_data_poly_true.t_bf_grid(1:end-1);...
    flip(G2_data_poly_true.t_bf_grid(2:end) - G2_data_poly_true.t_bf_grid(1:end-1))];

err1 = sum(abs(phi_1_poly_true - phi_1_poly_10dof).*weights1) / ...
    sum(abs(phi_1_poly_true).*weights1)

err2 = sum(abs(phi_2_poly_true - phi_2_poly_10dof).*weights2) / ...
    sum(abs(phi_2_poly_true).*weights2)


