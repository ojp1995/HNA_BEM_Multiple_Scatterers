% Test 3 HF iterative solve partial shadowing

clear all
addpath('../General_functions/')
addpath('../../BEAM_HNABEMLAB/')
addPathsHNA  % allows HNABEM to find all of the relevatn subfolders

vertices1 = [-2*pi 2*pi;
    0, 0];

vertices2 = [-4*pi -2*pi;
    2*pi 0];

R_max = 20;

kwave=10;
theta = 0;

Lgrad_coeff = 0.2;
alpha = 2;

C_wl_quad_outer = 1/10;

C_wl_quad_inner = 1/15;

[G1_data, G2_data, phi1_r, phi2_r, v_N1cell, v_N2cell, Xstruct1, Xstruct2] = ...
    HF_it_outer_function(kwave, vertices1, vertices2, R_max, theta, ...
    C_wl_quad_outer, C_wl_quad_inner, Lgrad_coeff, alpha);;

%% err calc

err1 = HF_it_L1_err_wrt_R(phi1_r{end}, phi1_r, R_max, G1_data.w_comb_outer);
err2 = HF_it_L1_err_wrt_R(phi2_r{end}, phi2_r, R_max, G2_data.w_comb_outer);

figure()
semilogy([0:2:(2*R_max - 3)], err1, 'DisplayName', '$L^{1}$ error of $\phi_{1}^{(r)}$')
hold on
semilogy([1:2:(2*R_max - 2)], err2, 'DisplayName', '$L^{1}$ error of $\phi_{2}^{(r)}$')
xlabel('r')
ylabel('Relative $L^{1}$ error')
legend show

%% plotting on bndy
figure()
for r = [1, 2, R_max]
    txt1 = ['r = ', mat2str(2*r-2)];
    plot(G1_data.t_mid_q_comb_outer/G1_data.L, real(phi1_r{r}), ...
        'DisplayName', txt1)
    hold on

end
xlim([-0.05 1.05])
ylim([-30 30])
xlabel('$x/L_{1}$')
ylabel('$\phi_{1}^{(r)}$')
title('HF iterative method $\phi_{1}^{(r)}$')
legend show

figure()
for r = [1, 2, R_max]
    txt2 = ['r = ', mat2str(2*r-1)];
    plot(G2_data.t_mid_q_comb_outer/G2_data.L, real(phi2_r{r}), ...
        'DisplayName', txt2)
    hold on

end
xlim([-0.05 1.05])
ylim([-30 30])
xlabel('$x/L_{2}$')
ylabel('$\phi_{2}^{(r)}$')
title('HF iterative method $\phi_{2}^{(r)}$')
legend show

%% plotting in domain

[u, ui, us] = HF_itproduce_plot_in_D(kwave, theta, G1_data, G2_data,...
    phi1_r{end}, phi2_r{end});