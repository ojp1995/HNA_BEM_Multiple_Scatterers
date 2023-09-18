function [err_normG1_it, err_normG2_it, aj_1_r, aj_2_r, aj1_1step, aj2_1step] = it_conv_test_fixed_dof(G1, G2, k, C1, C2, theta, C_wl, R_it, R_true) 
% In this function we will compute the error for a fixed degrees of freedom
% of increasing iterations approximation to the true solution with a high
% order number if iterations.
%
% Will be taking the l1 error.
%
% Problem parameters:
% (x1, y1) coordinates of \Gamma_{1}
% t1 paramterisation of \Gamma_{1}
% (x2, y2) coordinates of \Gamma_{2}
% t2 paramterisation of \Gamma_{2}
% k is the wavenumber
% C1, C2, constants for smoothing function
% theta, angle between downwards vertical and incident anti clockwise
% L1, L2, lengths of screens respecitively
%
% Discretisation paramters:
% C_wl is the number of dof per wavelength
% h1 is a vector of the step sizes for \Gamma_{1}
% h2 vector of step sizes for \Gamma_{2}

% N1 = ceil(k*L1./(C_wl*2*pi)) % number of itervals on G1
% N2 = ceil(k*L2./(C_wl*2*pi)) % number of intervals on G2

% computing variables for G1:
[x1, y1, t1, h1, h1vector, N1, L1] = discretisation_variables(G1, C_wl, k);
% computing variables for G2:
[x2, y2, t2, h2, h2vector, N2, L2] = discretisation_variables(G2, C_wl, k);

Q = 1;

[A1, x_col1 A1_sing, A1_smooth] = Rob_hop_screen_mat_poly_approx_PIM( x1, y1, x1, y1, h1vector, k, C1, C2, Q, t1);
[A2, x_col2, A2_sing, A2_smooth] = Rob_hop_screen_mat_poly_approx_PIM( x2, y2, x2, y2, h2vector, k, C1, C2, Q, t2);

S21 = S21_op_robhop(x1, y1, x2, y2, k, h1vector);
S12 = S21_op_robhop(x2, y2, x1, y1, k, h2vector);
% Solve in 1 step:
% constructing the matrix
A = [A1 S12 ;
    S21 A2];
% constructing the RHS
u_i_1step = zeros(length(h1vector)+length(h2vector), 1);
u_i_1step(1:length(h1vector), 1) = robhop_PW_incident(k, theta, x1, y1);
u_i_1step(length(h1vector)+1:length(h1vector)+length(h2vector), 1) = robhop_PW_incident(k, theta, x2, y2);
% doing the solve
aj_1step = A\u_i_1step;

% splitting up the coefficients
aj1_1step = aj_1step(1:length(h1vector), 1);
aj2_1step = aj_1step(length(h1vector)+1:length(h1vector)+length(h2vector), 1);
clear A u_i_1step aj_1step
%splitting the 

u_i_1_r = zeros(length(h1vector), R_it);
aj_1_r = zeros(length(h1vector), R_it);
u_i_2_r = zeros(length(h2vector), R_it);
aj_2_r = zeros(length(h2vector), R_it);


% r = 1 case:
u_i_1_r(:, 1) = robhop_PW_incident(k, theta, x1, y1);
aj_1_r(:, 1) = A1\u_i_1_r(:, 1);


u_i_2_r(:, 1) =  robhop_PW_incident(k, theta, x2, y2).' - S21*coeff_2soln_midpoint(aj_1_r(:, 1), L1, N1-1, N1);
aj_2_r(:, 1) = A2\u_i_2_r(:, 1);

% figure(15);
% plot(((1:N1) - 0.5)/N1, real(aj_1_r(:, 1)), 'DisplayName', 'r = 0')
% hold on
% 
% 
% figure(16);
% plot(((1:N2) - 0.5)/N2, real(aj_2_r(:, 1)), 'DisplayName', 'r = 1')
% hold on
% 
% figure(17)
% subplot(2, 1, 1)
% plot(((1:N1) - 0.5)/N1, aj_1_r(:, 1), 'DisplayName', sprintf('r = %g', 0))
% hold on
% subplot(2, 1, 2)
% plot(((1:N2) - 0.5)/N2, aj_2_r(:, 1), 'DisplayName', sprintf('r = %g', 1))
% hold on

% looping over the number of iterations
for r = 2:R_it
%     compute the even solutions on \Gamma_{1}
%     first compute the incident:
    u_i_1_r(:, r) =  robhop_PW_incident(k, theta, x1, y1).' - S12*coeff_2soln_midpoint(aj_2_r(:, r-1), L2, N2-1, N2);
%     conpute the coefficients
    aj_1_r(:, r) = A1\u_i_1_r(:, r);
   
%     figure(15)
%     plot(((1:N1) - 0.5)/N1, aj_1_r(:, r), 'DisplayName', sprintf('r = %g', 2*r - 2))
    
%     compute the odd solutions on \Gamma_{2}
%         first compute the incident:
    u_i_2_r(:, r) =  robhop_PW_incident(k, theta, x2, y2).' - S21*coeff_2soln_midpoint(aj_1_r(:, r), L1, N1-1, N1);
%     conpute the coefficients
    aj_2_r(:, r) = A2\u_i_2_r(:, r);
    
%      figure(16)
%     plot(((1:N2) - 0.5)/N2, aj_2_r(:, r), 'DisplayName', sprintf('r = %g', 2*r-1))
%     
%     figure(17)
%     subplot(2, 1, 1)
%     plot(((1:N1) - 0.5)/N1, aj_1_r(:, r), 'DisplayName', sprintf('r = %g', 2*r - 2))
%     
%     subplot(2, 1, 2)
%     plot(((1:N2) - 0.5)/N2, aj_2_r(:, r), 'DisplayName', sprintf('r = %g', 2*r-1))
%     
end


% figure(15)
% legend show
% title('Solution on \Gamma_{1} for increasing number of iterations')
% xlabel('x')
% ylabel('\phi_{1}^{r}')
% xlim([-0.1 1.1])
% figure(16)
% title('Solution on \Gamma_{2} for increasing number of iterations')
% xlabel('x')
% ylabel('\phi_{2}^{r}')
% xlim([-0.1 1.1])
% legend show
% 
% figure(17)
% subplot(2, 1, 1)
% % title('Solution on \Gamma_{1} for increasing number of iterations')
% xlabel('$s/L_{1}$', 'FontSize', 17, 'interpreter', 'latex')
% ylabel('Re $\left(\phi_{1}^{r} \right)$', 'FontSize', 17, 'interpreter', 'latex')
% xlim([-0.1 1.1])
% legend show
% set(legend,'fontsize',17);
% subplot(2, 1, 2)
% % title('Solution on \Gamma_{2} for increasing number of iterations')
% xlabel('$s/L_{2}$', 'FontSize', 17, 'interpreter', 'latex')
% ylabel('Re $ \left(\phi_{2}^{r} \right) $', 'FontSize', 17, 'interpreter', 'latex')
% xlim([-0.1 1.1])
% legend show
% set(legend,'fontsize',17);

    %%% computing higher order iterations
for r = R_it+1:R_true

    u_i_1_r(:, r) =  robhop_PW_incident(k, theta, x1, y1).' - S12*coeff_2soln_midpoint(aj_2_r(:, r-1), L2, N2-1, N2);
%     conpute the coefficients
    aj_1_r(:, r) = A1\u_i_1_r(:, r);

    %     compute the odd solutions on \Gamma_{2}
%         first compute the incident:
    u_i_2_r(:, r) =  robhop_PW_incident(k, theta, x2, y2).' - S21*coeff_2soln_midpoint(aj_1_r(:, r), L1, N1-1, N1);
%     conpute the coefficients
    aj_2_r(:, r) = A2\u_i_2_r(:, r);


     disp(r)
end


for r = 1:R_it
    % now computing the errors:
    err_normG1_it(r) = norm(aj_1_r(:, end) - aj_1_r(:, r), 1)/norm(aj_1_r(:, end), 1);
    err_normG2_it(r) = norm(aj_2_r(:, end) - aj_2_r(:, r), 1)/norm(aj_2_r(:, end), 1);

end
    
end
    