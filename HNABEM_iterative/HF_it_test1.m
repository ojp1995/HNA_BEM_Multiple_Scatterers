% Test 1 HF iterative solve

clear all
addpath('../General_functions/')
addpath('../../HNABEMLAB/')
addPathsHNA  % allows HNABEM to find all of the relevatn subfolders

vertices1 = [-2*pi 2*pi;
    0, 0];

vertices2 = [2*pi 0;
    5*pi 3*pi];

R_max = 20;

kwave= 10;
theta = 0;

C_wl_quad_outer = 1/10;

C_wl_quad_inner = 1/15;

Lgrad_coeff = 0.15;
alpha = 2;

tic
[G1_data, G2_data, phi1_r, phi2_r, v_N1cell, v_N2cell, Xstruct1, Xstruct2] = ...
    HF_it_outer_function_precompute_matrix(kwave, vertices1, vertices2, R_max, theta, ...
    C_wl_quad_outer, C_wl_quad_inner, Lgrad_coeff, alpha);
toc
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
title('$L^{1}$ error of $\phi_{j}^{(r)}$ with respect to the number of iterations, using the HNA iterative method - direct illumination')

%% plotting on bndy
figure()
for r =[1, 2, R_max]
    txt1 = ['r = ', mat2str(2*r-2)];
    plot(G1_data.t_mid_q_comb_outer/G1_data.L, real(phi1_r{r}), ...
        'DisplayName', txt1)
    hold on

end
xlim([-0.05 1.05])
ylim([-40 40])
xlabel('$x/L_{1}$')
ylabel('$\phi_{1}^{(r)}$')
title('HNA iterative method $\phi_{1}^{(r)}$ - direct illumination')
legend show

figure()
for r = [1, 2, R_max]
    txt2 = ['r = ', mat2str(2*r-1)];
    plot(G2_data.t_mid_q_comb_outer/G2_data.L, real(phi2_r{r}), ...
        'DisplayName', txt2)
    hold on

end
xlim([-0.05 1.05])
ylim([-40 40])
xlabel('$x/L_{2}$')
ylabel('$\phi_{2}^{(r)}$')
title('HNA iterative method $\phi_{2}^{(r)}$ - direct illumination')
legend show

%% plotting in domain
% 
% theta_rot = theta;
% 
% rot_mat = [cos(theta_rot) -sin(theta_rot);
%     sin(theta_rot) cos(theta_rot)];
% 
% 
% vert1x = rot_mat*[G1_data.G(1) G1_data.G(2)].';
% vert1y = rot_mat*[G1_data.G(3) G1_data.G(4)].';
% vertices1 = [vert1x.'; vert1y.'];
% 
% 
% vert2x = rot_mat*[G2_data.G(1) G2_data.G(2)].';
% vert2y = rot_mat*[G2_data.G(3) G2_data.G(4)].';
% vertices2 = [vert2x.'; vert2y.'];
% 
% G1_data.G = [vert1x.' vert1y.'];
% 
% G2_data.G = [vert2x.' vert2y.'];
% 
% G1_data_domain = get_graded_quad_points_HF_it(G1_data, C_wl_quad_outer,...
%     C_wl_quad_inner, kwave, Lgrad_coeff, alpha);
% 
% G2_data_domain = get_graded_quad_points_HF_it(G2_data, C_wl_quad_outer,...
%     C_wl_quad_inner, kwave, Lgrad_coeff, alpha);


[u, ui, us] = HF_itproduce_plot_in_D(kwave, theta, G1_data, G2_data,...
    phi1_r{end}, phi2_r{end});