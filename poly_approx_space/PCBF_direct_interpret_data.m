% load and interpret data from PCbf direct solver

clear all

addpath('../General_functions/')

load('PC_direct_test1_k10_theta0.mat')

G1_data.G = info_needed.G1;
G2_data.G = info_needed.G2;

k = info_needed.k;
theta = info_needed.theta;
Lgrad_coeff = info_needed.L_grad_coeff;
alpha = info_needed.alpha;
bf_dof_per_wl = info_needed.bf_dof_per_wl;

clear info_needed

G1_data = get_bf_graded_grid(G1_data, bf_dof_per_wl(1), k, ...
    Lgrad_coeff, alpha);

G2_data = get_bf_graded_grid(G2_data, bf_dof_per_wl(1), k, ...
    Lgrad_coeff, alpha);


t1_plot = linspace(0.05, G1_data.L/2, 3000);
t2_plot = linspace(0.05, G2_data.L/2, 3000);

for n = 1:length(bf_dof_per_wl)
    % first compute the mesh
    G1_data.G = [-2*pi, 2*pi, 0, 0];

    G2_data.G = [2*pi, 0, 5*pi, 3*pi]; 

    G1_data = get_bf_graded_grid(G1_data, bf_dof_per_wl(n), k, ...
        Lgrad_coeff, alpha);

    G2_data = get_bf_graded_grid(G2_data, bf_dof_per_wl(n), k, ...
        Lgrad_coeff, alpha);

    phi1(n, :) =  graded_coeff_2_solution(aj1_coeff{n}, ...
        G1_data.t_bf_grid, t1_plot, ...
        G1_data.L);

    phi2(n, :) =  graded_coeff_2_solution(aj2_coeff{n}, ...
        G2_data.t_bf_grid, t2_plot, ...
        G2_data.L);

end

t1_plot = [t1_plot.'; (G1_data.L - flip(t1_plot)).'];
t2_plot = [t2_plot.'; (G2_data.L - flip(t2_plot)).'];

% plot
figure()
for n = 1:length(bf_dof_per_wl)
    txt = ['dof = ', mat2str(1/bf_dof_per_wl(n)), 'per wavelength'];
    plot(t1_plot/G1_data.L, real(phi1(n, :)), 'DisplayName', txt)
    hold on

end
xlabel('$x/L_{1}$')
ylabel('$\phi_{1}(x)$')
title('Piecewise constant approximation to $\phi_{1}$')
legend show

figure()
for n = 1:length(bf_dof_per_wl)
    txt = ['dof = ', mat2str(1/bf_dof_per_wl(n)), 'per wavelength'];
    plot(t2_plot/G2_data.L, real(phi2(n, :)), 'DisplayName', txt)
    hold on

end
xlabel('$x/L_{2}$')
ylabel('$\phi_{2}(x)$')
title('Piecewise constant approximation to $\phi_{2}$')
legend show

%% quick and dirty error
for n = 1:length(bf_dof_per_wl)-1

    err_1(n) = sum(abs(phi1(n, :) - phi1(end, :)))/sum(abs(phi1(end, :)))

    err_2(n) = sum(abs(phi2(n, :) - phi2(end, :)))/sum(abs(phi2(end, :)))

end

figure()
semilogy(1./bf_dof_per_wl(1:end - 1), err_1, 'DisplayName', ...
    '$\ell_{1}$ error for $\phi_{1}$')
hold on
semilogy(1./bf_dof_per_wl(1:end-1), err_2, 'DisplayName', ...
    '$\ell_{1}$ error for $\phi_{2}$')
xlabel('Degrees of freedom per wavelength')
ylabel('$\ell_{1}$ error normalised')
title('Convergence test for piecewise constant basis function direct solver')