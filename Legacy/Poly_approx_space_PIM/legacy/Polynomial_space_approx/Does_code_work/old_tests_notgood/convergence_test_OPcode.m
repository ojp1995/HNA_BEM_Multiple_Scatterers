% this script is testing convergence for my own code to make sure that it
% is converging to the right thing.

clear all

addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/') % code to be tested
addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/Simon_code/') % simons code, taken to be as working code

load('NEWOPC_phi12_itall_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_n1')

% Nsample = ceil(1.5*(max(length(aj_1), length(aj_2))));
% N_sample1 = Nsample;
% N_sample2 = Nsample;
% labelling true solutions:
%iterative solutions
aj_1_r_true_it = aj_1_r;
aj_2_r_true_it = aj_2_r;
% for r =1:R
%     aj_1_r_true(:, r) = coeff_2soln_midpoint(aj_1_r(:, r), L1, N_sample1, N1);
%     aj_2_r_true(:, r) = coeff_2soln_midpoint(aj_2_r(:, r), L2, N_sample2,N2);
% end
clear aj_1_r aj_2_r

% solutions all in 1
aj_1_true_1step = aj_1;
aj_2_true_1step = aj_2;

clear aj_1 aj_2

%replace after rerun 
% [phi1, x_sample1] = coeff_2soln_midpoint(aj_1, L1, N1-1, N1);
% [phi2, x_sample2] = coeff_2soln_midpoint(aj_2, L2, N2-1, N2);



% All other information should be given from the loaded data so we now need
% to start the loop of testing

clear C_wl

C_wl = [1/5, 1/10, 1/20, 1/40, 1/80, 1/160]; % h <= C_wl \lambda, 1/C_wl intervals per wl.

N_sample1 = length(aj_1_r_true_it(:, 1));
h_sample1 = L1/N_sample1;
x_sample1 = [h_sample1/2: h_sample1: L1 - h_sample1/2];

N_sample2 = length(aj_2_r_true_it(:, 1));
h_sample2 = L2/N_sample2;
x_sample2 = [h_sample2/2: h_sample2: L2 - h_sample2/2];



figure(1)


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

    %%% plotting iteraitve solution
    figure(2*j-1)
    plot(x_col1/L1, real(aj_1_r(:, r)), 'DisplayName', 'r = 0')
    hold on
    title(['Solution on $\Gamma_{1}$, approx with ', num2str(1/C_wl(j)), 'dof per wl'])
    plot(x_sample1/L1, real(aj_1_r_true_it(:, end)), 'DisplayName', 'quasi true solution')
    
    figure(2*j)
    plot(x_col2/L2, real(aj_2_r(:, r)), 'DisplayName', 'r = 2')
    hold on
     title(['Solution on $\Gamma_{2}$, approx with ', num2str(1/C_wl(j)), 'dof per wl'])
    plot(x_sample2/L2, real(aj_2_r_true_it(:, end)), 'DisplayName', 'quasi true solution')
    % computing errors:
% %     err_G1(j, 1) = norm(abs( phi1 - phi1_r(:, 1), 2);
% %     err_G2(j, 1) = norm(abs( phi2 - phi2_r(:, 1), 2);
% % 
% %     err_norm_G1(j, 1) = norm(abs( phi1 - phi1_r), 2)/norm(phi1, 2);
% %     err_norm_G2(j, 1) = norm(abs( phi2 - phi2_r), 2)/norm(phi2, 2);

    % mapping "true" solution to same grid as current approximation
    % iterative true solution
    phi1_true_it = interp1(x_sample1, aj_1_r_true_it(:, r), x_col1);
    phi2_true_it = interp1(x_sample2, aj_2_r_true_it(:, r), x_col2);
    
    %%% 1 step true solution
    phi1_true_1step = interp1(x_sample1,  aj_1_true_1step, x_col1);
    phi2_true_1step = interp1(x_sample2,  aj_2_true_1step, x_col2);
    
    err_norm_it_interpl1(j, 2*r-1) = norm(phi1_true_it - aj_1_r(:, r).', 1)/norm(phi1_true_it, 1);
    err_norm_it_interpl1(j, 2*r) = norm(phi2_true_it - aj_2_r(:, r).', 1)/norm(phi2_true_it, 1);
    
    err_norm_it_interpl2(j, 2*r-1) = norm(phi1_true_it - aj_1_r(:, r).', 2)/norm(phi1_true_it, 2);
    err_norm_it_interpl2(j, 2*r) = norm(phi2_true_it - aj_2_r(:, r).', 2)/norm(phi2_true_it, 2);
    
    err_norm_it_interplinf(j, 2*r-1) = norm(phi1_true_it - aj_1_r(:, r).', inf)/norm(phi1_true_it, inf);
    err_norm_it_interplinf(j, 2*r) = norm(phi2_true_it - aj_2_r(:, r).', inf)/norm(phi2_true_it, inf);
    
    err_norm_it_interp_OP(j, 2*r-1) = sum(abs(phi1_true_it - aj_1_r(:, r).'))./sum(abs(phi1_true_it)); 
    err_norm_it_interp_OP(j, 2*r) = sum(abs(phi2_true_it - aj_2_r(:, r).'))./sum(abs(phi2_true_it)); 

    
%     err_norm_it(j, r) = norm(abs(aj_1_r_true(:, r) - phi1_r(:, r)), 1)/norm(aj_1_r_true(:, r), 1);
%     err_norm_it(j, r) = norm(abs(aj_2_r_true(:, r) - phi2_r(:, r)), 1)/norm(aj_2_r_true(:, r), 1);

%     figure(5);
%     plot(x_col1/L1, real(aj_1_r(:, 1)), 'DisplayName', 'r = 0')
%     hold on
% 
% 
%     figure(6);
%     plot(x_col2/L2, real(aj_2_r(:, 1)), 'DisplayName', 'r = 1')
%     hold on
% 
%     figure(7)
%     subplot(2, 1, 1)
%     plot(x_col1/L1, aj_1_r(:, 1), 'DisplayName', sprintf('r = %g', 0))
%     hold on
%     subplot(2, 1, 2)
%     plot(x_col2/L2, aj_2_r(:, 1), 'DisplayName', sprintf('r = %g', 1))
%     hold on

    %single step solve
    A = [A1 S12; S21 A2];
    u_inc = zeros(length(N1+N2), 1);
    u_inc(1:N1, 1) = PW_incident(k, theta, G1, x_col1);
    u_inc(N1+1:N1+N2, 1) = PW_incident(k, theta, G2, x_col2);
    
    aj_12 = A\u_inc;
    
    aj_1 = aj_12(1:N1, 1);
    aj_2 = aj_12(N1+1:N1+N2, 1);
    
    
    
    % true solution = iterative high dof, test solution, one step solve.
    err_norm_G12_1stepl1(j, 2*r-1) = norm( phi1_true_it - aj_1.' , 1)/norm(phi1_true_it, 1);
    err_norm_G12_1stepl1(j, 2*r) = norm( phi2_true_it - aj_2.' , 1)/norm(phi2_true_it, 1);
        
    err_norm_G12_1stepl2(j, 2*r-1) = norm( phi1_true_it - aj_1.' , 2)/norm(phi1_true_it, 2);
    err_norm_G12_1stepl2(j, 2*r) = norm( phi2_true_it - aj_2.' , 2)/norm(phi2_true_it, 2);
    
    % True solution one step method high dof, approx iterative method
    err_norm_G12_test_it_true_1stepl1(j, 2*r-1) = norm(phi1_true_1step - aj_1_r(:, r).', 1)/norm(phi1_true_1step, 1);
    err_norm_G12_test_it_true_1stepl1(j, 2*r) = norm(phi2_true_1step - aj_2_r(:, r).', 1)/norm(phi2_true_1step, 1);

    for r = 2:R
    %     compute the even solutions on \Gamma_{1}
    %     first compute the incident:
        u_i_1_r(:, r) =  PW_incident(k, theta, G1, x_col1).' - S12*coeff_2soln_midpoint(aj_2_r(:, r-1), L2, N2-1, N2);
    %     conpute the coefficients
        aj_1_r(:, r) = A1\u_i_1_r(:, r);
    %     compute the solution
%          [phi1_r(:, r), ~] = coeff_2soln_midpoint(aj_1_r(:, r), L1, N_sample1, N1);
%         figure(5)
%         plot(x_col1/L1, aj_1_r(:, r), 'DisplayName', sprintf('r = %g', 2*r - 2))

    %     compute the odd solutions on \Gamma_{2}
    %         first compute the incident:
        u_i_2_r(:, r) =  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_r(:, r), L1, N1-1, N1);
    %     conpute the coefficients
        aj_2_r(:, r) = A2\u_i_2_r(:, r);
    %     compute the solution
%          [phi2_r(:, r), ~] = coeff_2soln_midpoint(aj_2_r(:, r), L2, N_sample2, N2);

        % computing errors:
% %         err_G1(j, r) = norm(abs( phi1 - phi1_r(:, r) ), 2);
% %         err_G2(j, r) = norm(abs( phi2 - phi2_r(:, r) ), 2);
% % 
% %         err_norm_G1(j, r) = norm(abs( phi1 - phi1_r(:, r) ), 2)/norm(phi1, 2);
% %         err_norm_G2(j, r) = norm(abs( phi2 - phi2_r(:, r) ), 2)/norm(phi2, 2);
        phi1_true_it = interp1(x_sample1, aj_1_r_true_it(:, r), x_col1);
        phi2_true_it = interp1(x_sample2, aj_2_r_true_it(:, r), x_col2);

        err_norm_it_interpl1(j, 2*r-1) = norm(phi1_true_it - aj_1_r(:, r).', 1)/norm(phi1_true_it, 1);
        err_norm_it_interpl1(j, 2*r) = norm(phi2_true_it - aj_2_r(:, r).', 1)/norm(phi2_true_it, 1);
        
        err_norm_it_interpl2(j, 2*r-1) = norm(phi1_true_it - aj_1_r(:, r).', 2)/norm(phi1_true_it, 2);
        err_norm_it_interpl2(j, 2*r) = norm(phi2_true_it - aj_2_r(:, r).', 2)/norm(phi2_true_it, 2);
        
        err_norm_it_interplinf(j, 2*r-1) = norm(phi1_true_it - aj_1_r(:, r).', inf)/norm(phi1_true_it, inf);
        err_norm_it_interplinf(j, 2*r) = norm(phi2_true_it - aj_2_r(:, r).', inf)/norm(phi2_true_it, inf);
    
    
       err_norm_it_interp_OP(j, 2*r-1) = sum(abs(phi1_true_it - aj_1_r(:, r).'))./sum(abs(phi1_true_it)); 
       err_norm_it_interp_OP(j, 2*r) = sum(abs(phi2_true_it - aj_2_r(:, r).'))./sum(abs(phi2_true_it)); 

        
%         err_norm_it(j, r) = norm(abs(aj_1_r_true(:, r) - phi1_r(:, r)), 2)/norm(aj_1_r_true(:, r), 2);
%         err_norm_it(j, r) = norm(abs(aj_2_r_true(:, r) - phi2_r(:, r)), 2)/norm(aj_2_r_true(:, r), 2);
 
        % True solution one step method high dof, approx iterative method
        err_norm_G12_test_it_true_1stepl1(j, 2*r-1) = norm(phi1_true_1step - aj_1_r(:, r).', 1)/norm(phi1_true_1step, 1);
        err_norm_G12_test_it_true_1stepl1(j, 2*r) = norm(phi2_true_1step - aj_2_r(:, r).', 1)/norm(phi2_true_1step, 1);

%         
        figure(2*j-1)
        plot(x_col1/L1, real(aj_1_r(:, r)), 'DisplayName', sprintf('r = %g', 2*r - 2))
        legend show

        figure(2*j)
        plot(x_col2/L2, real(aj_2_r(:, r)), 'DisplayName', sprintf('r = %g', 2*r-1))
        legend show
        
        disp(r)
    end
    
end

err_norm_it_interpl1

err_norm_it_interp_OP

err_norm_it_interpl2

err_norm_it_interplinf

err_norm_G12_1stepl1

err_norm_G12_1stepl2

err_norm_G12_test_it_true_1stepl1

%%
%EOC computation for err_norm_itl1 qne G12 1 step l1
for r = 1:2*R
    for j = 1:(length(C_wl)-1)
        
        EOC_it_interpl1(j, r) = log2(err_norm_it_interpl1(j, r)/err_norm_it_interpl1(j+1, r));
        

    end
end

EOC_it_interpl1
% 
% for r = 1:2*R
%     for j = 1:(length(C_wl) - 1)
%         
%         EOC_allinone_l1(j, r) = log2(err_norm_G12_test_it_true_1stepl1(j, r)/err_norm_G12_test_it_true_1stepl1(j+1, r));
% 
%     end
%     
% end
%  
