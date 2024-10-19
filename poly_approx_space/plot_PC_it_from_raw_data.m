% loading data to then create plots

clear all


load('it_early_PC_direct_test5a_k10_thetapi_4.mat')

addpath('../General_functions/')

G1_data.G = info_needed.G1;
G2_data.G = info_needed.G2;

L_grad_coeff = info_needed.L_grad_coeff;
alpha = info_needed.alpha;

k = info_needed.k;

bf_dof_per_wl = info_needed.bf_dof_per_wl(1:end-1);


aj1_it_coeff = aj1_coeff;
aj2_it_coeff = aj2_coeff;

R_max = 10000;


G1_data = get_bf_graded_grid(G1_data, bf_dof_per_wl(end), k, ...
        L_grad_coeff, alpha);

G2_data = get_bf_graded_grid(G2_data, bf_dof_per_wl(end), k, ...
        L_grad_coeff, alpha);

x1_plotting = linspace(0.01, G1_data.L/2, 1000);
x2_plotting = linspace(0.01, G2_data.L/2, 1000);


C_wl_quad_err = 1/40;
G1_err_data = get_graded_quad_points(G1_data, C_wl_quad_err, k, ...
    L_grad_coeff, alpha);
G2_err_data = get_graded_quad_points(G2_data, C_wl_quad_err, k, ...
    L_grad_coeff, alpha);

% compute phi1 and phi2 at points of interest - for plotting
phi1_plotting = graded_coeff_2_solution(aj1_coeff{3}, ...
        G1_data.t_bf_grid, x1_plotting, ...
        G1_data.L);

phi2_plotting = graded_coeff_2_solution(aj2_coeff{3}, ...
        G2_data.t_bf_grid, x2_plotting, ...
        G2_data.L);

phi1_for_err = graded_coeff_2_solution(aj1_coeff{3}, ...
        G1_err_data.t_bf_grid, G1_err_data.t_mid_col, ...
        G1_data.L);

phi2_for_err = graded_coeff_2_solution(aj2_coeff{3}, ...
        G2_err_data.t_bf_grid, G2_err_data.t_mid_col, ...
        G2_data.L);
%%
it_of_interest = [1, 2, R_max];
figure()
for r = 1:length(it_of_interest)
    txt = ['r = ', mat2str(2*it_of_interest(r) - 2)];
    plot(G1_err_data.t_mid_col(10:end-10)/G1_data.L,...
        real(phi1_for_err(10: end - 10, it_of_interest(r))),...
        'DisplayName', txt)
    hold on

end
legend show
xlabel('$x/L_{1}$')
ylabel('$\phi_{1}^{(r)}$')
title('Iterative approximation to $\phi_{1}$ ')
xlim([-0.05 1.05])

figure()
for r = 1:length(it_of_interest)
    txt = ['r = ', mat2str(2*it_of_interest(r) - 1)];
    plot(G2_err_data.t_mid_col(10:end-10)/G2_data.L, ...
        real(phi2_for_err(10: end - 10, it_of_interest(r))),...
        'DisplayName', txt)
    hold on

end
legend show
xlabel('$x/L_{2}$')
ylabel('$\phi_{2}^{(r)}$')
title('Iterative approximation to $\phi_{2}$ ')
xlim([-0.05 1.05])


err_L1_G1 = L1_err_wrt_it_poly_it_bndy(G1_data, G1_err_data, aj_1_R{3}, R_max);

err_L1_G2 = L1_err_wrt_it_poly_it_bndy(G2_data, G2_err_data, aj_2_R{3}, R_max);


R_phi1 = [0:2:2*(R_max - 2)];
R_phi2 = [1:2:2*(R_max - 1)];

figure()
semilogy(R_phi1, err_L1_G1, 'DisplayName', '$\phi_{1}^{(r)}$ error')
hold on
semilogy(R_phi2, err_L1_G2, 'DisplayName', '$\phi_{2}^{(r)}$ error')
legend show
xlabel('Number of iterations, r')
ylabel('$\Vert \phi_{j}^{(R)} - \phi_{j}^{(r)} \Vert_{L^{1}((0, L_{j}))} / \Vert \phi_{j}^{(R)} \Vert_{L^{1}((0, L_{j}))}$')
title('$L^{1}$ error for an increasing number of iterations of $\phi_{j}^{(r)}$')
