% In this script we will be computing the l1 and L1 error for the solution
% on the boundary for the iterative solution POLYNOMIAL SOLVER ONLY
%
% script to get everything in line and then will functionise

clear all

% need to get some data to play with

Test1_poly_it
% would load otherwise
%%
C_wl_quad_err = 1/40;

G1_err_data = get_graded_quad_points(G1_data, C_wl_quad_err, k, ...
    Lgrad_coeff, alpha);

% compute phi1 for each iteration at each of the error calc points
phi1 = zeros(R_max, 2*length(G1_err_data.x_1_q));
figure()
for r = 1:R_max

    phi1(r, :) = graded_coeff_2_solution(aj_1_R(:, r), ...
        G1_data.t_bf_grid, G1_err_data.t_mid_q, G1_data.L);

    plot([ G1_err_data.t_mid_q; (G1_data.L - flip( G1_err_data.t_mid_q))], ...
        phi1(r, :))
    hold on

end

% Now computing the error
err_phi1 = zeros(R_max-1, 1);
for r = 1:R_max - 1

    err_phi1(r) = sum(abs(phi1(end, :) - phi1(r, :))./abs(phi1(end, :)))

end

figure()
semilogy(err_phi1(1:end))

%% Now looking L1 error

% phi1_normalised = phi1(end, :)*[G1_err_data.w; flip(G1_err_data.w)]
L1_err_G1 = zeros(R_max - 1, 1);
for r = 1:R_max-1

    L1_err_G1(r) = (abs(phi1(end, :) - phi1(r, :))./abs(phi1(end, :)))...
        *[G1_err_data.w; flip(G1_err_data.w)];

end
R_phi1 = [0:2:2*R_max - 4];
figure()
semilogy(R_phi1, L1_err_G1)
xlabel('Number of iterations, r')
ylabel('$\Vert \phi_{j}^{(R)} - \phi_{j}^{(r)} \Vert_{L^{1}((0, L_{j}))} / \Vert \phi_{j}^{(R)} \Vert_{L^{1}((0, L_{j}))}$')
title('$L^{1}$ error for an increasing number of iterations of $\phi_{j}^{(r)}$')

L1_err_wrt_it_poly_it_bndy(G1_data, G1_err_data, aj_1_R, R_max)
