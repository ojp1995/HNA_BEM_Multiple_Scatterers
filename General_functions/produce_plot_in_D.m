function [u, ui, us] = produce_plot_in_D(k, theta, G1_data, G2_data,...
    aj_1_R, aj_2_R)

x_min = min([G1_data.G(1), G1_data.G(3), G2_data.G(1), G2_data.G(3)]);
x_max = max([G1_data.G(1), G1_data.G(3), G2_data.G(1), G2_data.G(3)]);

y_min = min([G1_data.G(2), G1_data.G(4), G2_data.G(2), G2_data.G(4)]);
y_max = max([G1_data.G(2), G1_data.G(4), G2_data.G(2), G2_data.G(4)]);

% xy_min = min(x_min, y_min);
% xy_max = max(x_max, y_max);
% 
% X1 = -5 + xy_min: 0.1 :5 + xy_max;
% 
% X2 = -5 + xy_min: 0.1 :5 + xy_max;

X1 = -5 + x_min: 0.1 :5 + x_max;

X2 = -5 + y_min: 0.1 :5 + y_max;


[XX, YY] = meshgrid(X1, X2);

ui = incident(k, theta, XX, YY);

% collocation nodes and phi evaluated at those nodes!
phi_1 = graded_coeff_2_solution(aj_1_R(:, end), G1_data.t_bf_grid,...
    G1_data.t_mid_q, G1_data.L);

phi_2 = graded_coeff_2_solution(aj_2_R(:, end), G2_data.t_bf_grid,...
    G2_data.t_mid_q, G2_data.L);

us = soln_in_D_2_slow(G1_data, phi_1, G2_data, phi_2, k, X1, X2);

u = ui - us;

figure();
pcolor(XX, YY, real(u)); shading interp
hold on
plot([G1_data.G(1),G1_data.G(3)],[ G1_data.G(2),G1_data.G(4)], 'LineWidth', 3)
plot([G2_data.G(1),G2_data.G(3)],[ G2_data.G(2),G2_data.G(4)], 'LineWidth', 3)
shading interp
axis equal
colorbar
