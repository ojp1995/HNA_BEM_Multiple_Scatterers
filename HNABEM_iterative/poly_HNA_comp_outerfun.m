% direct comparison between polynomial iterative method and HNA method

% paths needed for HNA solve

clear all
addpath('../General_functions/') 
addpath('../../BEAM_HNABEMLAB/')
addPathsHNA  % allows HNABEM to find all of the relevatn subfolders

% Geometry
vertices1 = [-2*pi 2*pi;
    0, 0];

vertices2 = [2*pi 0;
    5*pi 3*pi];

R_max = 20;

kwave=10;
theta = 0;

C_wl_quad_outer = 1/20;

C_wl_quad_inner = 1/30;

Lgrad_coeff = 0.15;
alpha = 2;

%% HNA solve
[G1_data, G2_data, phi1_r, phi2_r, v_N1cell, v_N2cell, Xstruct1, Xstruct2] = ...
    HF_it_outer_function(kwave, vertices1, vertices2, R_max, theta, ...
    C_wl_quad_outer, C_wl_quad_inner, Lgrad_coeff, alpha);

%% poly solve
% adding relevant paths
addpath('../poly_approx_space/')

C1 = 1;
C2 = pi;

C_wl_bf1 = 1/20;
C_wl_bf2 = 1/20;

C_wl_quad= 1/20;

G1_data_poly.G = G1_data.G;
G2_data_poly.G = G2_data.G;

[aj1, aj2, uinc, G1_data_poly, G2_data_poly] = ...
    compute_2_screen_direct_poly(G1_data_poly, G2_data_poly, kwave, ...
    theta, C1, C2,  C_wl_bf1, C_wl_bf2, C_wl_quad, Lgrad_coeff, alpha);

% [G1_data_poly, G2_data_poly, aj_1_R, aj_2_R, us, phi_1_r, phi_2_r] = ...
%     compute_iteratuve_poly_scattering_prob_2_screens(G1_data_poly, ...
%     G2_data_poly, kwave, Lgrad_coeff, alpha, C_wl_bf1, C_wl_bf2, ...
%     C_wl_quad, R_max, theta, C1, C2, false, false);

%% Computing poly phi at the same points as 
phi_1_poly = graded_coeff_2_solution(aj1, ...
        G1_data_poly.t_bf_grid, G1_data.t_mid_q_outer, ...
        G1_data_poly.L);

phi_2_poly= graded_coeff_2_solution(aj2, ...
    G2_data_poly.t_bf_grid, G2_data.t_mid_q_outer, ...
    G2_data_poly.L);


% for r = 1:R_max
%     phi_1_poly(:, r) = graded_coeff_2_solution(aj_1_R(:, r), ...
%         G1_data_poly.t_bf_grid, G1_data.t_mid_q_outer, ...
%         G1_data_poly.L);
% 
%     phi_2_poly(:, r) = graded_coeff_2_solution(aj_2_R(:, r), ...
%         G2_data_poly.t_bf_grid, G2_data.t_mid_q_outer, ...
%         G2_data_poly.L);
%     
%     
% 
% end

%% plot to see a visual inspection
it_of_interest = [1, 2, R_max];
figure()
for r = 1:length(it_of_interest)
%     txt = ['Polynomial r = ', mat2str(2*it_of_interest(r) - 2)];
%     plot(G1_data.t_mid_q_comb_outer(10: end - 10)/G1_data.L,...
%         real(phi_1_poly(10: end - 10, it_of_interest(r))),...
%         'DisplayName', txt)
    hold on

    txt1 = ['HNA r = ', mat2str(2*it_of_interest(r) - 2)];
    plot(G1_data.t_mid_q_comb_outer(10: end - 10)/G1_data.L,...
        real(phi1_r{it_of_interest(r)}(10: end - 10)),...
        '--', 'DisplayName', txt1)

end
txt = 'Polynomial';
plot(G1_data.t_mid_q_comb_outer(10: end - 10)/G1_data.L,...
        real(phi_1_poly(10: end - 10)),...
        'DisplayName', txt)
legend show
xlabel('$x/L_{1}$')
ylabel('$\phi_{1}^{(r)}$')
title('Iterative approximation to $\phi_{1}$ ')
xlim([-0.05 1.05])

figure()
for r = 1:length(it_of_interest)
%     txt = ['Polynomial r = ', mat2str(2*it_of_interest(r) - 1)];
%     plot(G2_data.t_mid_q_comb_outer(10: end - 10)/G2_data.L, ...
%         real(phi_2_poly(10: end - 10, it_of_interest(r))),...
%         'DisplayName', txt)
    hold on

    txt1 = ['HNA r = ', mat2str(2*it_of_interest(r) - 1)];
    plot(G2_data.t_mid_q_comb_outer(10: end - 10)/G2_data.L, ...
        real(phi2_r{it_of_interest(r)}(10: end - 10)),...
        '--', 'DisplayName', txt1)

end
txt = 'Polynomial';
plot(G2_data.t_mid_q_comb_outer(10: end - 10)/G2_data.L, ...
        real(phi_2_poly(10: end - 10)),...
        'DisplayName', txt)
legend show
xlabel('$x/L_{2}$')
ylabel('$\phi_{2}^{(r)}$')
title('Iterative approximation to $\phi_{2}$ ')
xlim([-0.05 1.05])

%% error plots - poly is true soltuions

% polynomial computation
% C_wl_quad_err = 1/40;
% err_quadG1 = get_graded_quad_points(G1_data, C_wl_quad_err, kwave, ...
%     Lgrad_coeff, alpha);
% err_quadG2 = get_graded_quad_points(G2_data, C_wl_quad_err, kwave, ...
%     Lgrad_coeff, alpha);
% 
% phi_1_poly_for_err = graded_coeff_2_solution(aj1, ...
%         G1_data_poly.t_bf_grid, err_quadG1.t_mid_q_outer, ...
%         G1_data_poly.L);
% 
% phi_2_poly_for_err = graded_coeff_2_solution(aj2, ...
%     G2_data_poly.t_bf_grid, err_quadG2.t_mid_q_outer, ...
%     G2_data_poly.L);

err_phi1 = zeros(R_max - 1, 1);
err_phi2 = zeros(R_max - 1, 1);


for r = 1:R_max - 1

    err_phi1(r) = sum(abs(phi_1_poly - phi1_r{r})./abs(phi_1_poly)...
        .*G1_data.w_comb_outer);

    err_phi2(r) = sum(abs(phi_2_poly - phi2_r{r})./abs(phi_2_poly)...
        .*G2_data.w_comb_outer);

end

figure()
semilogy([0:2:(2*R_max - 3)], err_phi1, 'DisplayName', 'HNA $L^{1}$ error of $\phi_{1}^{(r)}$')
hold on
% semilogy([0:2:(2*R_max - 3)], err1_poly_HNA_true, 'DisplayName', 'Poly $L^{1}$ error of $\phi_{1}^{(r)}$')
% semilogy([0:2:(2*R_max - 3)], err_L1_G1, 'DisplayName', 'poly comp poly $\phi_{1}^{(r)}$ error')

% Gamma_2
semilogy([1:2:(2*R_max - 2)], err_phi2, 'DisplayName', 'HNA $L^{1}$ error of $\phi_{2}^{(r)}$')
% semilogy([1:2:(2*R_max - 2)], err2_poly_HNA_true, 'DisplayName', 'poly $L^{1}$ error of $\phi_{2}^{(r)}$')
% semilogy([1:2:(2*R_max - 2)], err_L1_G2, 'DisplayName', 'poly comp poly $\phi_{2}^{(r)}$ error')

xlabel('r')
ylabel('Relative $L^{1}$ error, poly solver true solution')
legend show



%% error plots - taking HNA as true solution

err1_HNA = HF_it_L1_err_wrt_R(phi1_r{end}, phi1_r, R_max, ...
    G1_data.w_comb_outer);
err2_HNA = HF_it_L1_err_wrt_R(phi2_r{end}, phi2_r, R_max, ...
    G2_data.w_comb_outer);

% polynomial computation
C_wl_quad_err = 1/40;
err_quadG1 = get_graded_quad_points(G1_data, C_wl_quad_err, kwave, ...
    Lgrad_coeff, alpha);
err_quadG2 = get_graded_quad_points(G2_data, C_wl_quad_err, kwave, ...
    Lgrad_coeff, alpha);


err1_poly_HNA_true = zeros(R_max - 1, 1);
err2_poly_HNA_true = zeros(R_max - 1, 1);

err_L1_G1 = L1_err_wrt_it_poly_it_bndy(G1_data_poly, err_quadG1, aj_1_R, R_max);

err_L1_G2 = L1_err_wrt_it_poly_it_bndy(G2_data_poly, err_quadG2, aj_2_R, R_max);



for r = 1:R_max - 1

    err1_poly_HNA_true(r) = sum(abs( phi1_r{end}(10: end - 10) - phi_1_r(10: end - 10, r) )...
        ./abs(phi1_r{end}(10: end - 10)).*G1_data.w_comb_outer(10: end - 10));

    err2_poly_HNA_true(r) = sum(abs( phi2_r{end}(10: end  - 10) - phi_2_r(10: end  - 10, r) )...
        ./abs(phi2_r{end}(10: end  - 10)).*G2_data.w_comb_outer(10: end  - 10));

end


figure()
semilogy([0:2:(2*R_max - 3)], err1_HNA, 'DisplayName', 'HNA $L^{1}$ error of $\phi_{1}^{(r)}$')
hold on
semilogy([0:2:(2*R_max - 3)], err1_poly_HNA_true, 'DisplayName', 'Poly $L^{1}$ error of $\phi_{1}^{(r)}$')
semilogy([0:2:(2*R_max - 3)], err_L1_G1, 'DisplayName', 'poly comp poly $\phi_{1}^{(r)}$ error')

% Gamma_2
semilogy([1:2:(2*R_max - 2)], err2_HNA, 'DisplayName', 'HNA $L^{1}$ error of $\phi_{2}^{(r)}$')
semilogy([1:2:(2*R_max - 2)], err2_poly_HNA_true, 'DisplayName', 'poly $L^{1}$ error of $\phi_{2}^{(r)}$')
semilogy([1:2:(2*R_max - 2)], err_L1_G2, 'DisplayName', 'poly comp poly $\phi_{2}^{(r)}$ error')

xlabel('r')
ylabel('Relative $L^{1}$ error')
legend show





