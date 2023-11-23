function [G1_data, G2_data, aj_1_R, aj_2_R, us] = ...
    compute_iteratuve_poly_scattering_prob_2_screens(G1_data, G2_data, ...
    k, Lgrad_coeff, alpha, C_wl_bf1, C_wl_bf2, C_wl_quad, R_max, theta,...
    C1, C2, bndy_plot, Domain_plot)
% In this function we compute the scattering problem using the iterative
% polynomial solver using a graded midpoint and graded product midpoint
% rule. 
%
% Input parameters:
% [G1_Data, G2_Data], initially the coordinates of the screen 
% Gamma1 and Gamma 2
%
% k, wavenumber
% [Lgrad_coeff, alpha], parameters for grading of grids
% C_wl_bf1, C_wl_bf2, num of dof per wavelength for basis fucntions
% C_wl_quad, number of dof for quadrature
% R_max, number of iterations
% Optional (ish):
% bndy_plot, Domain_plot, if either true, plots will be given and soln in
% domain computed.
% 
%
% Output:
% [G1_data, G2_data], all information relating to screens, bf_grid,
% quad_grid, collocation points
% [aj_1_R, aj_2_R], coefficients
% us, if Domain_plot = true, then scattered field output, otherwise 0 value
% given as output.


% computing the basis functions
G1_data = get_bf_graded_grid(G1_data, C_wl_bf1, k, Lgrad_coeff, alpha);

G2_data = get_bf_graded_grid(G2_data, C_wl_bf2, k, Lgrad_coeff, ...
    alpha);

% quadrature nodes and other info needed
G1_data = get_graded_quad_points(G1_data, C_wl_quad, k, ...
    Lgrad_coeff, alpha);

G2_data = get_graded_quad_points(G2_data, C_wl_quad, k, ...
    Lgrad_coeff, alpha);

% collocation points
G1_data = manipulate_collocation_points_graded(G1_data);
G2_data = manipulate_collocation_points_graded(G2_data);

% computing the matrix
[S11, S12, S21, S22, u_inc1, u_inc2] = ...
    compute_matrices_for_iterative_solve(G1_data, G2_data, k, ...
    G1_data.t_bf_grid, G2_data.t_bf_grid, theta, C1, C2 );

% iterative solve
[aj_1_R, aj_2_R, phi_1_r, phi_2_r] = iterative_poly_graded_PIM_solve(...
    S11, S12, S21, S22, u_inc1, u_inc2, R_max, G1_data.t_bf_grid, ...
    G2_data.t_bf_grid, G1_data, G2_data);

% produce plots of soln on bndy - Optional argument

if bndy_plot == true
    figure();
    for r = 1:R_max
        txt1 = ['r = ', mat2str(2*r-2)];
        plot(G1_data.s/G1_data.L, phi_1_r(:, r), 'DisplayName', txt1);
        hold on
    
    end
    legend show
    xlabel('$x/L_{1}$')
    ylabel('$\phi_{1}^{(r)}$')
    title('Iterative approximation to $\phi_{1}$ ')
    
    figure();
    for r = 1:R_max
        txt1 = ['r = ', mat2str(2*r - 1)];
        plot(G2_data.s/G2_data.L, phi_2_r(:, r), 'DisplayName', txt1);
        hold on
    
    end
    legend show
    xlabel('$x/L_{2}$')
    ylabel('$\phi_{2}^{(r)}$')
    title('Iterative approximation to $\phi_{2}$ ')
end

% produce plot in Domain - Optional argument

if Domain_plot == true
    [~, ~, us] = produce_plot_in_D(k ,theta, G1_data, G2_data, aj_1_R,...
    aj_2_R);

else
    us = 0;

end



