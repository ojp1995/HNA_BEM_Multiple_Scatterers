% outer function to compute convergence of two screens using the piecewise
% constant basis functions direct solver, i.e., not an iterative method

% Test 2 - iterative

clear all

addpath('../General_functions/')  % access to solvers needed

% geometry set up
G1_data.G = [0, 0, 3*pi, 2*pi];

G2_data.G = [pi, -2, 3*pi, -1]; 

% coefficients needed for creating grid for basis functions and quadrature
% points
Lgrad_coeff = 0.15;
alpha = 2;

k = 10;

theta = 3*pi/8;
C1 = 1;
C2 = pi;

bf_dof_per_wl = [1/5, 1/10, 1/20, 1/40, 1/80];

% creating empty cells for coefficients
aj1_coeff = {};
aj2_coeff = {};


for n = 1:length(bf_dof_per_wl)

    G1_data.G = [0, 0, 3*pi, 2*pi];

    G2_data.G = [pi, -2, 3*pi, -1];  
    
    tic
    [G1_data, G2_data, aj1_coeff{n}, aj2_coeff{n}, ~, ~, ~] = ...
    compute_iteratuve_poly_scattering_prob_2_screens(G1_data, G2_data, ...
    k, Lgrad_coeff, alpha, bf_dof_per_wl(n), bf_dof_per_wl(n), ...
    bf_dof_per_wl(n), R_max, theta, C1, C2, false, false);

    toc

    clear G1_data G2_data

end

%% Plotting and computing error

% Create mesh we want to evaluate points at

G1_data.G = [0, 0, 3*pi, 2*pi];

G2_data.G = [pi, -2, 3*pi, -1]; 

G1_data = get_bf_graded_grid(G1_data, bf_dof_per_wl(1), k, ...
    Lgrad_coeff, alpha);

G2_data = get_bf_graded_grid(G2_data, bf_dof_per_wl(1), k, ...
    Lgrad_coeff, alpha);


t1_plot = linspace(0.05, G1_data.L/2, 3000);
t2_plot = linspace(0.05, G2_data.L/2, 3000);

for n = 1:length(bf_dof_per_wl)
    % first compute the mesh
    G1_data.G = [0, 0, 3*pi, 2*pi];

    G2_data.G = [pi, -2, 3*pi, -1];  

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

%%
% quick save only the bare minimum
info_needed.G1 = G1_data.G;
info_needed.G2 = G2_data.G;
info_needed.k = k;
info_needed.theta = theta;
info_needed.L_grad_coeff = Lgrad_coeff;
info_needed.alpha = alpha;
info_needed.bf_dof_per_wl = bf_dof_per_wl;

save('it_PC_direct_test2_k10_theta3pi_8', 'aj1_coeff', 'aj2_coeff', 'info_needed')






