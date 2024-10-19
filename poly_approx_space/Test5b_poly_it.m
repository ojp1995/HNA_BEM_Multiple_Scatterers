% test 5a - near trapping, almost touching screens

clear all

addpath('../General_functions/')  % access to solvers needed

% introducing the screens, storing the data in a struct object 

G1_data.G = [0, 0, pi, 0];

G2_data.G = [1/2, 3/2, pi, 1]; 

% coefficients needed for creating grid for basis functions and quadrature
% points
Lgrad_coeff = 0.15;
alpha = 2;

% creating basis function information
C_wl_bf1 = 1/10;
C_wl_bf2 = 1/10;

C_wl_quad= 1/20;

R_max = 150;

k = 10;

theta = pi/4;
C1 = 1;
C2 = pi;


% solve
[G1_data, G2_data, aj_1_R, aj_2_R, us, phi_1_r, phi_2_r] = ...
    compute_iteratuve_poly_scattering_prob_2_screens(G1_data, G2_data, ...
    k, Lgrad_coeff, alpha, C_wl_bf1, C_wl_bf2, C_wl_quad, R_max, theta, ...
    C1, C2, false, true);

%% Plot of relevant iterative solutions on bndy
it_of_interest = [1, 2, R_max];
figure()
for r = 1:length(it_of_interest)
    txt = ['r = ', mat2str(2*it_of_interest(r) - 2)];
    plot(G1_data.s(10:end-10)/G1_data.L,...
        real(phi_1_r(10: end - 10, it_of_interest(r))),...
        'DisplayName', txt)
    hold on

end
legend show
xlabel('$x/L_{1}$')
ylabel('$\phi_{1}^{(r)}$')
title('PC iterative approximation to $\phi_{1}$ - funnel screens')
xlim([-0.05 1.05])

figure()
for r = 1:length(it_of_interest)
    txt = ['r = ', mat2str(2*it_of_interest(r) - 1)];
    plot(G2_data.s(10:end-10)/G2_data.L, ...
        real(phi_2_r(10: end - 10, it_of_interest(r))),...
        'DisplayName', txt)
    hold on

end
legend show
xlabel('$x/L_{2}$')
ylabel('$\phi_{2}^{(r)}$')
title('PC iterative approximation to $\phi_{1}$ - funnel screens')
xlim([-0.05 1.05])

%% compute error wrt to iterations
C_wl_quad_err = 1/40;
G1_err_data = get_graded_quad_points(G1_data, C_wl_quad_err, k, ...
    Lgrad_coeff, alpha);
G2_err_data = get_graded_quad_points(G2_data, C_wl_quad_err, k, ...
    Lgrad_coeff, alpha);

err_L1_G1 = L1_err_wrt_it_poly_it_bndy(G1_data, G1_err_data, aj_1_R, R_max);

err_L1_G2 = L1_err_wrt_it_poly_it_bndy(G2_data, G2_err_data, aj_2_R, R_max);

R_phi1 = [0:2:2*(R_max - 2)];
R_phi2 = [1:2:2*(R_max - 1)];

figure()
semilogy(R_phi1, err_L1_G1, 'DisplayName', '$\phi_{1}^{(r)}$ error')
hold on
semilogy(R_phi2, err_L1_G2, 'DisplayName', '$\phi_{2}^{(r)}$ error')
legend show
xlabel('Number of iterations, r')
ylabel('$\Vert \phi_{j}^{(R)} - \phi_{j}^{(r)} \Vert_{L^{1}((0, L_{j}))} / \Vert \phi_{j}^{(R)} \Vert_{L^{1}((0, L_{j}))}$')
title('$L^{1}$ error of $\phi_{j}^{(r)}$ with respect to the number of iterations, using the PC iterative method - funnel screens')
