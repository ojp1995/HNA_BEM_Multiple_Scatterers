% looking at the difference between the single screen solve and the
% iterative method

% loading high dof single screen solve

clear all

figure(1)
figure(2)
close figure 1
close figure 2

addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/') % code to be tested
addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/Simon_code/') % simons code, taken to be as working code

load('NEWOPC_phi12_itall_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_n1')

%clearing iterative solutions
clear aj_1_r aj_2_r

% solutions all in 1
aj_1_true_1step = aj_1;
aj_2_true_1step = aj_2;

clear aj_1 aj_2


clear C_wl

C_wl = 1/20;

R = 5

N_sample1 = length(aj_1_true_1step(:, 1));
h_sample1 = L1/N_sample1;
x_sample1 = [h_sample1/2: h_sample1: L1 - h_sample1/2];

N_sample2 = length(aj_2_true_1step(:, 1));
h_sample2 = L2/N_sample2;
x_sample2 = [h_sample2/2: h_sample2: L2 - h_sample2/2];

for j = 1:length(C_wl)
    r = 1;
    % number of intervals on wach screen:
     N1 = ceil(k*L1./(C_wl(j)*2*pi)) % number of itervals on G1
     N2 = ceil(k*L2./(C_wl(j)*2*pi)) % number of intervals on G2



    Q = 1;  % quadrature discretisation
    % iterative solve
    
    [A1, x_col1 A1_sing, A1_smooth] = screen_mat_poly_approx_PIM(N1, 0, L1, k, C1, C2, Q);
    [A2, x_col2, A2_sing, A2_smooth] = screen_mat_poly_approx_PIM(N2, 0, L2, k, C1, C2, Q);

    [S21, col_points_1] = S21_op(G1, G2, k, N1, N2);
    [S12, col_points_2] = S21_op(G2, G1, k, N2, N1);
    % automating the re-reflections

    % initialising
    u_i_1_r = zeros(length(x_col1), R);
    aj_1_r = zeros(length(x_col1), R);
    u_i_2_r = zeros(length(x_col2), R);
    aj_2_r = zeros(length(x_col2), R);
%     phi1_r = zeros(N_sample1+1, R);
%     phi2_r = zeros(N_sample2+1, R);
    % r = 1 case:
    u_i_1_r(:, 1) = PW_incident(k, theta, G1, x_col1);
    aj_1_r(:, 1) = A1\u_i_1_r(:, 1);
%     [phi1_r(:, 1), x_sample1] = coeff_2soln_midpoint(aj_1_r(:, 1), L1, N_sample1, N1);
    

    u_i_2_r(:, 1) =  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_r(:, 1), L1, N1-1, N1);
    aj_2_r(:, 1) = A2\u_i_2_r(:, 1);
%     [phi2_r(:, 1), x_sample2] = coeff_2soln_midpoint(aj_2_r(:, 1), L2, N_sample2, N2);

      %%% 1 step true solution
    phi1_true_1step = interp1(x_sample1,  aj_1_true_1step, x_col1);
    phi2_true_1step = interp1(x_sample2,  aj_2_true_1step, x_col2);
    
    figure(2*j-1)
    plot(x_col1/L1, real(aj_1_r(:, r)), 'DisplayName', 'r = 0')
    hold on
    title(['Solution on $\Gamma_{1}$, approx with ', num2str(1/C_wl(j)), 'dof per wl'])
    
    figure(2*j)
    plot(x_col2/L2, real(aj_2_r(:, r)), 'DisplayName', 'r = 1')
    hold on
    title(['Solution on $\Gamma_{2}$, approx with ', num2str(1/C_wl(j)), 'dof per wl'])
    
    % plotting difference between approximation and truth
    approx_diff_1 = aj_1_r(:, r) - phi1_true_1step.';
    figure(4*R + 2*j-1)
    plot(x_col1, real(approx_diff_1), 'DisplayName', 'r = 0')
    hold on
    title(['Difference in approx on $\Gamma_{1}$, approx with ', num2str(1/C_wl(j)), 'dof per wl'])
    
    approx_diff_2 = aj_2_r(:, r) - phi2_true_1step.';
    figure(4*R + 2*j)
    plot(x_col2, real(approx_diff_2), 'DisplayName', 'r = 1')
    hold on
    title(['Difference in approx on $\Gamma_{2}$, approx with ', num2str(1/C_wl(j)), 'dof per wl'])
    
    
     err_norm_G12_test_it_true_1stepl1(j, 2*r-1) = norm(phi1_true_1step - aj_1_r(:, r).', 1)/norm(phi1_true_1step, 1);
     err_norm_G12_test_it_true_1stepl1(j, 2*r) = norm(phi2_true_1step - aj_2_r(:, r).', 1)/norm(phi2_true_1step, 1);

     for r = 2:R
    %     compute the even solutions on \Gamma_{1}
    %     first compute the incident:
        u_i_1_r(:, r) =  PW_incident(k, theta, G1, x_col1).' - S12*coeff_2soln_midpoint(aj_2_r(:, r-1), L2, N2-1, N2);
    %     conpute the coefficients
        aj_1_r(:, r) = A1\u_i_1_r(:, r);
        figure(1)
        plot(x_col1/L1, real(aj_1_r(:, r)), 'DisplayName', sprintf('r = %g', 2*r - 2))
        legend show
        
    %     
    %     compute the odd solutions on \Gamma_{2}
    %         first compute the incident:
        u_i_2_r(:, r) =  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_r(:, r), L1, N1-1, N1);
    %     conpute the coefficients
        aj_2_r(:, r) = A2\u_i_2_r(:, r);
        figure(2)
         plot(x_col2/L2, real(aj_2_r(:, r)), 'DisplayName', sprintf('r = %g', 2*r-1))
        legend show
        
        approx_diff_1 = aj_1_r(:, r) - phi1_true_1step.';
        
        figure(4*R + 2*j-1)
        plot(x_col1, real(approx_diff_1), 'DisplayName', sprintf('r = %g', 2*r - 2))
        legend show
        
        approx_diff_2 = aj_2_r(:, r) - phi2_true_1step.';
        
        figure(4*R + 2*j)
        plot(x_col2, real(approx_diff_2), 'DisplayName', sprintf('r = %g', 2*r-1))
        legend show
        
        err_norm_G12_test_it_true_1stepl1(j, 2*r-1) = norm(phi1_true_1step - aj_1_r(:, r).', 1)/norm(phi1_true_1step, 1);
     err_norm_G12_test_it_true_1stepl1(j, 2*r) = norm(phi2_true_1step - aj_2_r(:, r).', 1)/norm(phi2_true_1step, 1);
    
    
     end
     
     figure(2*j-1)
     plot(x_sample1/L1, real(aj_1_true_1step), 'DisplayName', 'All in 1 solution')
     
     figure(2*j)
     plot(x_sample2/L2, real(aj_2_true_1step), 'DisplayName', 'All in 1 solution')
end
