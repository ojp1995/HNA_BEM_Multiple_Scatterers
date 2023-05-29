% in this code we will be looking at the convergence of the iterative
% method with respect to the number of iterations. we are taking the true
% solution to be 99th iteration for the first screen and the 100th
% iteration for the second screen. (It might be slightly shifted I don't
% know, I know I let R = 50.) 

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This method does not quite work, the dof dominates and you can't get any
% closer to the true solution after about 2/3 iterations. In an new
% function it will be tested to see if for a fixed number of dofs that the
% iterations converge.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% adding paths and loading data
addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/Does_code_work/old_tests_notgood') % code to be tested
addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx')

load('NEWOPC_phi12_itall_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_n1')

% all in 1 solution
phi1_1step = aj_1;
phi2_1step = aj_2;
clear aj_1 aj_2


% relabelling the true solutions
%iterative solutions
aj_1_r_true_it = aj_1_r;
aj_2_r_true_it = aj_2_r;

clear aj_1_r aj_2_r

% % solutions all in 1
% aj_1_true_1step = aj_1;
% aj_2_true_1step = aj_2;
% 
% clear aj_1 aj_2

clear C_wl
%%
C_wl = [1/5, 1/10, 1/20, 1/40, 1/80, 1/160] %, 1/320, 1/640]; % h <= C_wl \lambda, 1/C_wl intervals per wl.

N_sample1 = length(aj_1_r_true_it(:, 1));
h_sample1 = L1/N_sample1;
x_sample1 = [h_sample1/2: h_sample1: L1 - h_sample1/2];

N_sample2 = length(aj_2_r_true_it(:, 1));
h_sample2 = L2/N_sample2;
x_sample2 = [h_sample2/2: h_sample2: L2 - h_sample2/2];

clear R
R = 10;
for j = 1:length(C_wl)
    r = 1;
    disp(r)
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
    
    % solutions for \Gamma_{1}
    u_i_1_r(:, 1) = PW_incident(k, theta, G1, x_col1);
    aj_1_r(:, 1) = A1\u_i_1_r(:, 1);
    
    % solutions on \Gamma_{2}
    u_i_2_r(:, 1) =  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_r(:, 1), L1, N1-1, N1);
    aj_2_r(:, 1) = A2\u_i_2_r(:, 1);
    
    % computing errors
    
    % mapping "true" solution to same grid as current approximation
    % iterative true solution
    phi1_true_it_final = interp1(x_sample1, aj_1_r_true_it(:, end), x_col1);
    phi2_true_it_final = interp1(x_sample2, aj_2_r_true_it(:, end), x_col2);
    
    % all in one solution
    phi1_all1_true = interp1(x_sample1, phi1_1step, x_col1);
    phi2_all1_true = interp1(x_sample2, phi2_1step, x_col2);
    
    %%% plotting iteraitve solution
    figure(2*j-1)
    plot(x_col1/L1, real(aj_1_r(:, r)), 'DisplayName', 'r = 0')
    hold on
    title(['Solution on $\Gamma_{1}$, approx with ', num2str(1/C_wl(j)), 'dof per wl'])
    plot(x_col1/L1, real(phi1_true_it_final), ':*', 'DisplayName', 'quasi true solution interpolate')
%     plot(x_sample1/L1, real(aj_1_r_true_it(:, end)), ':*', 'DisplayName', 'quasi true solution')
    
    figure(2*j)
    plot(x_col2/L2, real(aj_2_r(:, r)), 'DisplayName', 'r = 2')
    hold on
    title(['Solution on $\Gamma_{2}$, approx with ', num2str(1/C_wl(j)), 'dof per wl'])
     plot(x_col2/L2, real(phi2_true_it_final), ':*', 'DisplayName', 'quasi true solution interpolate')
%     plot(x_sample2/L2, real(aj_2_r_true_it(:, end)),':*',  'DisplayName', 'quasi true solution')

    
    
    
    err_normG1_it_interpl1(r, j) = norm(phi1_true_it_final - aj_1_r(:, r).', 1)/norm(phi1_true_it_final, 1);
    err_normG2_it_interpl1(r, j) = norm(phi2_true_it_final - aj_2_r(:, r).', 1)/norm(phi2_true_it_final, 1);
        
%     err_normG1_it_interpl1_awayends(j, r) = norm(phi1_true_it_final(2:end-1) - aj_1_r(2:end-1, r).', 1)/norm(phi1_true_it_final(2:end-1), 1);
%     err_normG2_it_interpl1_awayends(j, r) = norm(phi2_true_it_final(2:end-1) - aj_2_r(2:end-1, r).', 1)/norm(phi2_true_it_final(2:end-1), 1);
%     
    err_normG1_1step_comp(r, j) = norm(phi1_all1_true - aj_1_r(:, r).', 1)/norm(phi1_all1_true, 1);
    err_normG2_1step_comp(r, j) = norm(phi2_all1_true - aj_2_r(:, r).', 1)/norm(phi2_all1_true, 1);
    % now begin loop over iterations
    for r = 2:R
    %     compute the even solutions on \Gamma_{1}
    %     first compute the incident:
        u_i_1_r(:, r) =  PW_incident(k, theta, G1, x_col1).' - S12*coeff_2soln_midpoint(aj_2_r(:, r-1), L2, N2-1, N2);
    %     conpute the coefficients
        aj_1_r(:, r) = A1\u_i_1_r(:, r);
        
        %     compute the odd solutions on \Gamma_{2}
    %         first compute the incident:
        u_i_2_r(:, r) =  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_r(:, r), L1, N1-1, N1);
    %     conpute the coefficients
        aj_2_r(:, r) = A2\u_i_2_r(:, r);
        
        err_normG1_it_interpl1(r, j) = norm(phi1_true_it_final - aj_1_r(:, r).', 1)/norm(phi1_true_it_final, 1);
        err_normG2_it_interpl1(r, j) = norm(phi2_true_it_final - aj_2_r(:, r).', 1)/norm(phi2_true_it_final, 1);
        
%         err_normG1_it_interpl1_awayends(j, r) = norm(phi1_true_it_final(2:end-1) - aj_1_r(2:end-1, r).', 1)/norm(phi1_true_it_final(2:end-1), 1);
%         err_normG2_it_interpl1_awayends(j, r) = norm(phi2_true_it_final(2:end-1) - aj_2_r(2:end-1, r).', 1)/norm(phi2_true_it_final(2:end-1), 1);
    
        err_normG1_1step_comp(r, j) = norm(phi1_all1_true - aj_1_r(:, r).', 1)/norm(phi1_all1_true, 1);
        err_normG2_1step_comp(r, j) = norm(phi2_all1_true - aj_2_r(:, r).', 1)/norm(phi2_all1_true, 1);
    
        figure(2*j-1)
        plot(x_col1/L1, real(aj_1_r(:, r)), 'DisplayName', sprintf('r = %g', 2*r - 2))
        legend show

        figure(2*j)
        plot(x_col2/L2, real(aj_2_r(:, r)), 'DisplayName', sprintf('r = %g', 2*r-1))
        legend show
        
        disp(r)
        
    end
end

err_normG1_it_interpl1

err_normG2_it_interpl1

err_normG1_1step_comp

err_normG2_1step_comp

%%
% computing the EOC for each screen w.r.t r

for j = 1:length(C_wl)-1
    
    for r = 1:R
        
%         EOC_G1_wrt_r(r, j) = log2(err_normG1_it_interpl1(r, j)/err_normG1_it_interpl1(r+1, j));
%         
%         EOC_G2_wrt_r(r, j) = log2(err_normG2_it_interpl1(r, j)/err_normG2_it_interpl1(r+1, j));
%         
        
        EOC_G1_wrt_r(j, r) = log2(err_normG1_it_interpl1(r, j)/err_normG1_it_interpl1(r, j+1));
        
        EOC_G2_wrt_r(j, r) = log2(err_normG2_it_interpl1(r, j)/err_normG2_it_interpl1(r, j+1));
        
        EOC_G1_allin1 (j, r) = log2( err_normG1_1step_comp(r, j)/err_normG1_1step_comp(r, j+1) );
        
        EOC_G2_allin1 (j, r) = log2( err_normG2_1step_comp(r, j)/err_normG2_1step_comp(r, j+1) );

    end
end

EOC_G1_wrt_r

EOC_G2_wrt_r

%% plotting the errors
x_G1_err_plots = [0:2:2*R-1];
figure()
hold on
for j = 1:length(C_wl)
    
    semilogy(x_G1_err_plots, err_normG1_it_interpl1(:, j), 'LineStyle', '-.', 'DisplayName', sprintf(' dof = %g', 1/C_wl(j)) )
    
    title('Approximation of $\phi_{1}^{(r)}$ on $\Gamma_{1}$ to the solution with $N=640$ degrees of freedom where r = 50 (change)')
end
xlabel('r')
ylabel('Normalised $\ell_{1}$ error')

legend show

x_G2_err_plots = [1:2:2*R];
figure()
hold on
for j = 1:length(C_wl)
    
    semilogy(x_G2_err_plots, err_normG2_it_interpl1(:, j), 'LineStyle', '-.', 'DisplayName', sprintf(' dof = %g', 1/C_wl(j)) )
    
    title('Approximation of $\phi_{2}^{(r)}$ on $\Gamma_{2}$ to the solution with $N=640$ degrees of freedom where r = 50 (change)')
end
xlabel('r')
ylabel('Normalised $\ell_{1}$ error')

legend show



