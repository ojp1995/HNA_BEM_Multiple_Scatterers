% Test 1 - In this test we will be looking at two screens that directly
% illuminate each other

clear all

addpath('../General_functions/')  % access to solvers needed
tic
% introducing the screens, storing the data in a struct object 

G1_data.G = [-2*pi, 2*pi, 0, 0];

G2_data.G = [2*pi, 0, 5*pi, 3*pi]; 

% coefficients needed for creating grid for basis functions and quadrature
% points
Lgrad_coeff = 0.15;
alpha = 1;

% creating basis function information
C_wl_bf1 = 1/2;
C_wl_bf2 = 1/2;

C_wl_quad= 1/20;

R_max = 20;

k = 10;

theta = 0;
C1 = 1;
C2 = pi;



% solve
[G1_data, G2_data, aj_1_R, aj_2_R, us, phi_1_r, phi_2_r] = ...
    compute_iteratuve_poly_scattering_prob_2_screens(G1_data, G2_data, ...
    k, Lgrad_coeff, alpha, C_wl_bf1, C_wl_bf2, C_wl_quad, R_max, theta, ...
    C1, C2, false, false);

toc

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
title('PC iterative approximation to $\phi_{1}$ - direct illumination')
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
title('PC iterative approximation to $\phi_{2}$ - direct illumination')
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
title('$L^{1}$ error of $\phi_{j}^{(r)}$ with respect to the number of iterations, using the PC iterative method - direct illumination')

% [u, ui, us] = produce_plot_in_D(k, theta, G1_data, G2_data,...
%     aj_1_R(:, end), aj_2_R(:, end));


% G1_data = get_bf_graded_grid(G1_data, C_wl_bf1, k, Lgrad_coeff, alpha);
% 
% G2_data = get_bf_graded_grid(G2_data, C_wl_bf2, k, Lgrad_coeff, ...
%     alpha);
% 
% % creating quadrature points
% 
% % number of degrees of freedom per wavelength for the quadrature
% 
% % quadrature nodes and other information needed
% G1_data = get_graded_quad_points(G1_data, C_wl_quad, k, ...
%     Lgrad_coeff, alpha);
% 
% G2_data = get_graded_quad_points(G2_data, C_wl_quad, k, ...
%     Lgrad_coeff, alpha);
% 
% % collocation points
% G1_data = manipulate_collocation_points_graded(G1_data);
% G2_data = manipulate_collocation_points_graded(G2_data);
% 
% [S11, S12, S21, S22, u_inc1, u_inc2] = ...
%     compute_matrices_for_iterative_solve(G1_data, G2_data, k, ...
%     G1_data.t_bf_grid, G2_data.t_bf_grid, theta, C1, C2 );
% 
% % iterative solve
% 
% 
% [aj_1_R, aj_2_R, phi_1_r, phi_2_r] = iterative_poly_graded_PIM_solve(...
%     S11, S12, S21, S22, u_inc1, u_inc2, R_max, G1_data.t_bf_grid, ...
%     G2_data.t_bf_grid, G1_data, G2_data);
% 
% % produce plots of soln on bndy - Optional argument
% 
% if bndy_plot == true
%     figure();
%     for r = 1:R_max
%         txt1 = ['r = ', mat2str(2*r-2)];
%         plot(G1_data.s/G1_data.L, phi_1_r(:, r), 'DisplayName', txt1);
%         hold on
%     
%     end
%     legend show
%     xlabel('$x/L_{1}$')
%     ylabel('$\phi_{1}^{(r)}$')
%     title('Iterative approximation to $\phi_{1}$ ')
%     
%     figure();
%     for r = 1:R_max
%         txt1 = ['r = ', mat2str(2*r - 1)];
%         plot(G2_data.s/G2_data.L, phi_2_r(:, r), 'DisplayName', txt1);
%         hold on
%     
%     end
%     legend show
%     xlabel('$x/L_{2}$')
%     ylabel('$\phi_{2}^{(r)}$')
%     title('Iterative approximation to $\phi_{2}$ ')
% end
% 
% % produce plot in Domain - Optional argument
% % 
% Domain_plot= true
% if Domain_plot == true
% [u, ui, us] = produce_plot_in_D(k, theta, G1_data, G2_data,...
%     aj_1_R(:, end), aj_2_R(:, end));

% end

