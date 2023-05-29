%parallel screen example

clear all

tic

% % introducing the screens
% len = 2;
% d = 2*pi; % distance between screens
% 
% G1 = [-len*pi, len*pi, 0, 0];
% 
% G2 = [d - len*pi, len*pi, d, 0];

C_wl= 1/40

k = 5;  % wavenumber

theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

% number of iterations for approximation and truth
R_it = 100;
R_true = 200;

len = [0.5, 1, 2, 3, 4, 5, 6, 8, 10];
d = 2*pi; % distance between screens

for j = 1: length(len)
    
    % introducing the screens


    G1 = [-len(j)*pi, len(j)*pi, 0, 0];

    G2 = [d - len(j)*pi, len(j)*pi, d, 0];
    
    [err_normG1_it(:, j), err_normG2_it(:, j), ~, ~] = it_conv_test_fixed_dof( G1, G2, k, C1, C2, theta, C_wl, R_it, R_true); 


    
end


%%
x_G1_err_plots = [0:2:2*R_it-1];
figure()
hold on
for j = 1:length(len)
    G1 = [-len(j)*pi, len(j)*pi, 0, 0];
    L1 = sqrt( (G1(3) - G1(1))^2 + (G1(4) - G1(2))^2 ); 
    semilogy(x_G1_err_plots, err_normG1_it(:, j), 'LineStyle', '-.', 'DisplayName', sprintf(' L1 = %g',L1) )
end
title('Approximation of $\phi_{1}^{(r)}$ on $\Gamma_{1}$ to the "true solution" with the same degrees of freedom as the approximation where r is much larger')
set(gca, 'YScale', 'log') % But you can explicitly force it to be logarithmic
xlabel('r')
ylabel('Normalised $\ell_{1}$ error')
xlim([-0.1 x_G1_err_plots(end)+0.1])
% ylim([-0.1 0.8])
legend show

%%
x_G2_err_plots = [1:2:2*R_it];
figure()
hold on
for j = 1:length(len)
    G2 = [d - len(j)*pi, len(j)*pi, d, 0];
    L2 = sqrt( (G2(3) - G2(1))^2 + (G2(4) - G2(2))^2 ); 
    semilogy(x_G2_err_plots, err_normG2_it(:, j), 'LineStyle', '-.', 'DisplayName', sprintf(' L2 = %g', L2) )
end
title('Approximation of $\phi_{2}^{(r)}$ on $\Gamma_{2}$ to the "true solution" with the same degrees of freedom as the approximation where r is much larger ')
set(gca, 'YScale', 'log') % But you can explicitly force it to be logarithmic
% grid on
xlabel('r')
ylabel('Normalised $\ell_{1}$ error')
xlim([-0.1 x_G2_err_plots(end)+0.1])
% ylim([-0.0001 0.02])
legend show

% % % % N1 = ceil(k*L1./(C_wl*2*pi)) % number of itervals on G1
% % % % N2 = ceil(k*L2./(C_wl*2*pi)) % number of intervals on G2
% % % % 
% % % % % step size
% % % % h1 = L1/N1;
% % % % h2 = L2/N2;
% % % % 
% % % % % computing the midpoints coordinates for \Gamma_{1} and \Gamma_{2}
% % % % x1 = G1(1)+((1:N1)-0.5)*(G1(3)-G1(1))/N1;
% % % % y1 = G1(2)+((1:N1)-0.5)*(G1(4)-G1(2))/N1;
% % % % 
% % % % % what if I added a paramterisation here for the screen?
% % % % t1 = [0:h1:L1];
% % % % 
% % % % x2 = G2(1)+((1:N2)-0.5)*(G2(3)-G2(1))/N2;
% % % % y2 = G2(2)+((1:N2)-0.5)*(G2(4)-G2(2))/N2;
% % % % 
% % % % % what if I added a paramterisation here for the screen?
% % % % t2 = [0:h2:L2];
% % % % 
% % % % h1vector = h1*ones(size(x1)); % Lengths of the elements on screen 1 
% % % % h2vector = h2*ones(size(x2)); % and on screen 2
% % % % h = [h1vector,h2vector]; % the lengths of the elements
% % % % x = [x1,x2]; % the x-coords of the element midpoints
% % % % y = [y1,y2]; % and their y-coords

% % % Q = 1;  % quadrature discretisation
% % % warning('v1_mid_weights.m has been bodged so it can work, need to look more closely')
% % % warning('Line 41 of robhop_midpoint_m2_phi.m, should it bw h or h_new???')
% % % [A1, x_col1 A1_sing, A1_smooth] = Rob_hop_screen_mat_poly_approx_PIM( x1, y1, x1, y1, h1vector, k, C1, C2, Q, t1);
% % % [A2, x_col2, A2_sing, A2_smooth] = Rob_hop_screen_mat_poly_approx_PIM( x2, y2, x2, y2, h2vector, k, C1, C2, Q, t2);
% % % 
% % % S21 = S21_op_robhop(x1, y1, x2, y2, k, h1vector);
% % % S12 = S21_op_robhop(x2, y2, x1, y1, k, h2vector);
% % % 
% % % %%
% % % % compute the solutions from coeff2soln
% % % N_sample = 1000;  % number of points we are sampling the solution at
% % % 
% % % 
% % % % automating the re-reflections
% % % R = 10; % this is the number of re-reflections we want
% % % % initialising
% % % u_i_1_r = zeros(length(h1vector), R);
% % % aj_1_r = zeros(length(h1vector), R);
% % % u_i_2_r = zeros(length(h2vector), R);
% % % aj_2_r = zeros(length(h2vector), R);
% % % phi1_r = zeros(N_sample, R);
% % % phi2_r = zeros(N_sample, R);
% % % % r = 1 case:
% % % u_i_1_r(:, 1) = robhop_PW_incident(k, theta, x1, y1);
% % % aj_1_r(:, 1) = A1\u_i_1_r(:, 1);
% % % [phi1_r(:, 1), x_sample1] = coeff_2soln_midpoint(aj_1_r(:, 1), L1, N_sample-1, N1);
% % % 
% % % figure();
% % % plot(((1:N1) - 0.5)/N1, real(aj_1_r(:, 1)), 'DisplayName', 'r = 0')
% % % hold on
% % % 
% % % u_i_2_r(:, 1) =  robhop_PW_incident(k, theta, x2, y2).' - S21*coeff_2soln_midpoint(aj_1_r(:, 1), L1, N1-1, N1);
% % % aj_2_r(:, 1) = A2\u_i_2_r(:, 1);
% % % [phi2_r(:, 1), x_sample2] = coeff_2soln_midpoint(aj_2_r(:, 1), L2, N_sample-1, N2);
% % % 
% % % figure(15);
% % % plot(((1:N1) - 0.5)/N1, real(aj_1_r(:, 1)), 'DisplayName', 'r = 0')
% % % hold on
% % % 
% % % 
% % % figure(16);
% % % plot(((1:N2) - 0.5)/N2, real(aj_2_r(:, 1)), 'DisplayName', 'r = 1')
% % % hold on
% % % 
% % % figure(17)
% % % subplot(2, 1, 1)
% % % plot(((1:N1) - 0.5)/N1, aj_1_r(:, 1), 'DisplayName', sprintf('r = %g', 0))
% % % hold on
% % % subplot(2, 1, 2)
% % % plot(((1:N2) - 0.5)/N2, aj_2_r(:, 1), 'DisplayName', sprintf('r = %g', 1))
% % % hold on
% % % % 
% % % for r = 2:R
% % % %     compute the even solutions on \Gamma_{1}
% % % %     first compute the incident:
% % %     u_i_1_r(:, r) =  robhop_PW_incident(k, theta, x1, y1).' - S12*coeff_2soln_midpoint(aj_2_r(:, r-1), L2, N2-1, N2);
% % % %     conpute the coefficients
% % %     aj_1_r(:, r) = A1\u_i_1_r(:, r);
% % % %     compute the solution
% % %     [phi1_r(:, r), ~] = coeff_2soln_midpoint(aj_1_r(:, r), L1, N_sample-1, N1);
% % %     figure(15)
% % %     plot(((1:N1) - 0.5)/N1, aj_1_r(:, r), 'DisplayName', sprintf('r = %g', 2*r - 2))
% % %     
% % % %     compute the odd solutions on \Gamma_{2}
% % % %         first compute the incident:
% % %     u_i_2_r(:, r) =  robhop_PW_incident(k, theta, x2, y2).' - S21*coeff_2soln_midpoint(aj_1_r(:, r), L1, N1-1, N1);
% % % %     conpute the coefficients
% % %     aj_2_r(:, r) = A2\u_i_2_r(:, r);
% % % %     compute the solution
% % %     [phi2_r(:, r), ~] = coeff_2soln_midpoint(aj_2_r(:, r), L2, N_sample-1, N2);
% % %     
% % %      figure(16)
% % %     plot(((1:N2) - 0.5)/N2, aj_2_r(:, r), 'DisplayName', sprintf('r = %g', 2*r-1))
% % %     
% % %     figure(17)
% % %     subplot(2, 1, 1)
% % %     plot(((1:N1) - 0.5)/N1, aj_1_r(:, r), 'DisplayName', sprintf('r = %g', 2*r - 2))
% % %     
% % %     subplot(2, 1, 2)
% % %     plot(((1:N2) - 0.5)/N2, aj_2_r(:, r), 'DisplayName', sprintf('r = %g', 2*r-1))
% % %     
% % % end
% % % 
% % % 
% % % figure(15)
% % % legend show
% % % title('Solution on \Gamma_{1} for increasing number of iterations')
% % % xlabel('x')
% % % ylabel('\phi_{1}^{r}')
% % % xlim([-0.1 1.1])
% % % figure(16)
% % % title('Solution on \Gamma_{2} for increasing number of iterations')
% % % xlabel('x')
% % % ylabel('\phi_{2}^{r}')
% % % xlim([-0.1 1.1])
% % % legend show
% % % 
% % % figure(17)
% % % subplot(2, 1, 1)
% % % % title('Solution on \Gamma_{1} for increasing number of iterations')
% % % xlabel('$s/L_{1}$', 'FontSize', 17, 'interpreter', 'latex')
% % % ylabel('Re $\left(\phi_{1}^{r} \right)$', 'FontSize', 17, 'interpreter', 'latex')
% % % xlim([-0.1 1.1])
% % % legend show
% % % set(legend,'fontsize',17);
% % % subplot(2, 1, 2)
% % % % title('Solution on \Gamma_{2} for increasing number of iterations')
% % % xlabel('$s/L_{2}$', 'FontSize', 17, 'interpreter', 'latex')
% % % ylabel('Re $ \left(\phi_{2}^{r} \right) $', 'FontSize', 17, 'interpreter', 'latex')
% % % xlim([-0.1 1.1])
% % % legend show
% % % set(legend,'fontsize',17);
