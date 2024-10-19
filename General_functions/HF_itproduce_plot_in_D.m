function [u, ui, us] = HF_itproduce_plot_in_D(k, theta, G1_data, G2_data,...
    phi1, phi2)
% Plot in the domain for HF iterative method

x_min = min([G1_data.G(1), G1_data.G(3), G2_data.G(1), G2_data.G(3)]);
x_max = max([G1_data.G(1), G1_data.G(3), G2_data.G(1), G2_data.G(3)]);

y_min = min([G1_data.G(2), G1_data.G(4), G2_data.G(2), G2_data.G(4)]);
y_max = max([G1_data.G(2), G1_data.G(4), G2_data.G(2), G2_data.G(4)]);

X1 = -5 + x_min: 0.1 :5 + x_max;

X2 = -5 + y_min: 0.1 :5 + y_max;

[XX, YY] = meshgrid(X1, X2);

ui = incident(k, theta, XX, YY);

us = HF_it_soln_in_2D(G1_data, phi1, G2_data, phi2, k, X1, X2);

u = ui - us;

figure();
pcolor(XX, YY, real(u)); shading interp
hold on
plot([G1_data.G(1),G1_data.G(3)],[ G1_data.G(2),G1_data.G(4)], 'LineWidth', 3)
plot([G2_data.G(1),G2_data.G(3)],[ G2_data.G(2),G2_data.G(4)], 'LineWidth', 3)
shading interp
axis equal
colorbar