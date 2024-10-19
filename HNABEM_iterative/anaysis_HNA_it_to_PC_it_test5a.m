% computing relative error between HNA iterative and PC iterative

clear all

addpath('../General_functions/')
addpath('..//poly_approx_space/')

% load in the HNA solution

load('test5a_L2pi_6kit_10_march.mat')

R_max = length(phi1_r);


% load poly solver
load('it_PC_direct_test5a_k10_thetapi_4.mat')

%%
% evaluate the polynomial solver at the points the HNA it is evaluated at
poly_select = 2;
G1_data_poly.G = info_needed.G1;
G2_data_poly.G = info_needed.G2;
Lgrad_coeff = info_needed.L_grad_coeff;
alpha = info_needed.alpha;
C_wl_poly = info_needed.bf_dof_per_wl(poly_select);
k = info_needed.k;


G1_data_poly = get_bf_graded_grid(G1_data_poly, C_wl_poly, k, Lgrad_coeff, alpha);

G2_data_poly = get_bf_graded_grid(G2_data_poly, C_wl_poly, k, Lgrad_coeff, ...
    alpha);

phi1_poly = graded_coeff_2_solution(aj1_coeff{end}(:, end), ...
    G1_data_poly.t_bf_grid, ...
    G1_data.t_mid_q_outer, G1_data.L);


phi2_poly = graded_coeff_2_solution(aj2_coeff{poly_select}(:, end), ...
    G2_data_poly.t_bf_grid, ...
    G2_data.t_mid_q_outer, G2_data.L);

for r = 1:R_max

    err1(r) = sum( ((abs(phi1_poly - phi1_r{r})./abs(phi1_poly)))...
        .*G1_data.w_comb_outer);

    err2(r) = sum( ((abs(phi2_poly - phi2_r{r})./abs(phi2_poly)))...
        .*G2_data.w_comb_outer);

end

r1 = [0:2:2*R_max-2];
r2 = [1:2: 2*R_max - 1];


%%
figure()
semilogy(r1, err1, 'DisplayName','$\phi_{1}$')
hold on
semilogy(r2, err2, 'DisplayName','$\phi_{2}$')
legend show

xlabel('r')
ylabel('Realative error')

title('Relative error between the HNA iterative method and the final iteration of the piecewise constant method for the case of parallel screens.')

xlim([-200, 2*R_max+ 200])

