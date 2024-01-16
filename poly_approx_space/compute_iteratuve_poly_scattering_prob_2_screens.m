function [G1_data, G2_data, aj_1_R, aj_2_R, us, phi_1_r, phi_2_r] = ...
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


% % information for aligning grids
% tol1 = max([G1_data.w; 1e-2]);
% tol2 = max([G2_data.w; 1e-2]);
% 
% 
% [G1grid_quad] = align_quad_grid(G1_data.t_bf_grid, G1_data.t_grid, tol1);
% 
% [G2grid_quad] = align_quad_grid(G2_data.t_bf_grid, G2_data.t_grid, tol2);
% 
% % Now need to get the nodes from this grid and the weights and the 2d
% % versions as well
% 
% % G1 midpoints (parameterised) and weights
% for j = 1:(length(G1grid_quad)-1)
%     G1t_mid(j) = (G1grid_quad(j+1) + G1grid_quad(j))/2;
%     G1w(j) = G1grid_quad(j+1) - G1grid_quad(j);
% end
% 
% % constructing 2D points
% G1x_1_1 = G1_data.G(1) + G1t_mid*(G1_data.G(3) - G1_data.G(1))/G1_data.L;
% G1y_1_1 = G1_data.G(2) + G1t_mid*(G1_data.G(4) - G1_data.G(2))/G1_data.L;
% 
% G1x_1_2 = G1_data.G(1) + (G1_data.L - G1t_mid)*(G1_data.G(3) - G1_data.G(1))/G1_data.L;
% G1y_1_2 = G1_data.G(2) + (G1_data.L - G1t_mid)*(G1_data.G(4) - G1_data.G(2))/G1_data.L;
% 
% % G2 midpoints (parameterised) and weights
% for j = 1:(length(G2grid_quad)-1)
%     G2t_mid(j) = (G2grid_quad(j+1) + G2grid_quad(j))/2;
%     G2w(j) = G2grid_quad(j+1) - G2grid_quad(j);
% end
% 
% % constructing 2D points
% G2x_1_1 = G2_data.G(1) + G2t_mid*(G2_data.G(3) - G2_data.G(1))/G2_data.L;
% G2y_1_1 = G2_data.G(2) + G2t_mid*(G2_data.G(4) - G2_data.G(2))/G2_data.L;
% 
% G2x_1_2 = G2_data.G(1) + (G2_data.L - G2t_mid)*(G2_data.G(3) - G2_data.G(1))/G2_data.L;
% G2y_1_2 = G2_data.G(2) + (G2_data.L - G2t_mid)*(G2_data.G(4) - G2_data.G(2))/G2_data.L;
% 
% 
% G1_data.t_grid = G1grid_quad;
% G1_data.t_mid_q = G1t_mid;
% G1_data.w = G1w;
% G1_data.x_1_q = G1x_1_1;
% G1_data.y_1_q = G1y_1_1;
% G1_data.x_2_q = G1x_1_2;
% G1_data.y_2_q = G1y_1_2;
% 
% 
% G2_data.t_grid = G2grid_quad;
% G2_data.t_mid_q = G2t_mid;
% G2_data.w = G2w;
% G2_data.x_1_q = G2x_1_1;
% G2_data.y_1_q = G2y_1_1;
% G2_data.x_2_q = G2x_1_2;
% G2_data.y_2_q = G2y_1_2;


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
        plot(G1_data.s(10:end-10)/G1_data.L, phi_1_r(10:end-10, r), 'DisplayName', txt1);
        hold on
    
    end
    legend show
    xlabel('$x/L_{1}$')
    ylabel('$\phi_{1}^{(r)}$')
    title('Iterative approximation to $\phi_{1}$ ')
    xlim([-0.05 1.05])
    
    figure();
    for r = 1:R_max
        txt1 = ['r = ', mat2str(2*r - 1)];
        plot(G2_data.s(10:end-10)/G2_data.L, phi_2_r(10:end-10, r), 'DisplayName', txt1);
        hold on
        
    
    end
    legend show
    xlabel('$x/L_{2}$')
    ylabel('$\phi_{2}^{(r)}$')
    title('Iterative approximation to $\phi_{2}$ ')
    xlim([-0.05 1.05])
end

% produce plot in Domain - Optional argument

if Domain_plot == true
    [~, ~, us] = produce_plot_in_D(k ,theta, G1_data, G2_data, aj_1_R,...
    aj_2_R);

else
    us = 0;

end



