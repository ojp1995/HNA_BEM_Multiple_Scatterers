% in this script we are testing to see if for a fixed number of degrees of
% freedom the iterative solver converges to a higher order iteration.

clear all

% adding paths and loading data
addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/Does_code_work/old_tests_notgood') % code to be tested
addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx')

G1 = [-2*pi, 2*pi, 2*pi, -2*pi];

G2 = [pi, 2*pi, 5*pi/2, pi]; 
% G1 = [ -2*pi, 2*pi,0, 0]; 
% G1 = [-pi/2, 2*pi, 0, 0];  % tighter formation
L1 = sqrt( (G1(3) - G1(1))^2 + (G1(4) - G1(2))^2 );  % length of G1
% G2 = [4*pi, 0, 7*pi, 3*pi];
% G2 = [pi/2, 0, pi, 3*pi];
L2 = sqrt( (G2(3) - G2(1))^2 + (G2(4) - G2(2))^2 );  % length of G2

k = 5;  % wavenumber

theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

% I think the set up has now been introduced so we can start to compute the
% matrices for \Gamma_{1} and \Gamma_{2}

Q = 1;  % quadrature discretisation

C_wl = [1/5, 1/10, 1/20, 1/40, 1/80, 1/160, 1/320]; % h <= C_wl \lambda, 1/C_wl intervals per wl.

R_it = 15;

R_true = 30;  % this is the number of iterations in the true case.


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
    
    
    
    % string condition numbers (crazy ideas but could be worth while)
    cond_A_G1(j) = cond(A1);
    
    cond_A_G2(j) = cond(A2);
    
    cond_A_S21(j) = cond(S21);
    cond_A_S12(j) = cond(S12);

    % initialising
    u_i_1_r = zeros(length(x_col1), R_true);
    aj_1_r = zeros(length(x_col1), R_true);
    u_i_2_r = zeros(length(x_col2), R_true);
    aj_2_r = zeros(length(x_col2), R_true);
    
     % solutions for \Gamma_{1}
    u_i_1_r(:, 1) = PW_incident(k, theta, G1, x_col1);
    aj_1_r(:, 1) = A1\u_i_1_r(:, 1);
    
    % solutions on \Gamma_{2}
    u_i_2_r(:, 1) =  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_r(:, 1), L1, N1-1, N1);
    aj_2_r(:, 1) = A2\u_i_2_r(:, 1);
    
    figure(2*j-1)
    plot(x_col1/L1, real(aj_1_r(:, r)), 'DisplayName', 'r = 0')
    hold on
    
    figure(2*j)
    plot(x_col2/L2, real(aj_2_r(:, r)), 'DisplayName', 'r = 1')
    hold on
    
    % now begin loop over iterations
    for r = 2:R_it
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
        
        
        figure(2*j-1)
        plot(x_col1/L1, real(aj_1_r(:, r)), 'DisplayName', sprintf('r = %g', 2*r - 2))
        

        figure(2*j)
        plot(x_col2/L2, real(aj_2_r(:, r)), 'DisplayName', sprintf('r = %g', 2*r-1))
        
        
        disp(r)
        
    end
    %%% computing higher order iterations
    for r = R_it+1:R_true
        
        u_i_1_r(:, r) =  PW_incident(k, theta, G1, x_col1).' - S12*coeff_2soln_midpoint(aj_2_r(:, r-1), L2, N2-1, N2);
    %     conpute the coefficients
        aj_1_r(:, r) = A1\u_i_1_r(:, r);
        
        %     compute the odd solutions on \Gamma_{2}
    %         first compute the incident:
        u_i_2_r(:, r) =  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_r(:, r), L1, N1-1, N1);
    %     conpute the coefficients
        aj_2_r(:, r) = A2\u_i_2_r(:, r);
        

         disp(r)
    end
    
    figure(100)
    plot(x_col1/L1, real(aj_1_r(:, end)), 'DisplayName', sprintf('N = %g', 1/C_wl(j)))
    xlabel('x/L')
    ylabel('$\phi_{1}^{(98)}$')
    title('Comparison of different orders of degrees of freedom approximation to 98th iteration on $\Gamma_{1}$')
    hold on
    
    figure(101)
    plot(x_col2/L2, real(aj_2_r(:, end)), 'DisplayName', sprintf('N = %g', 1/C_wl(j)))
    xlabel('x/L')
    ylabel('$\phi_{1}^{(99)}$')
    title('Comparison of different orders of degrees of freedom approximation to 99th iteration on $\Gamma_{2}$')
    hold on
    
    for r = 1:R_it
    % now computing the errors:
        err_normG1_it_interpl1(r, j) = norm(aj_1_r(:, end) - aj_1_r(:, r), 1)/norm(aj_1_r(:, end), 1);
        err_normG2_it_interpl1(r, j) = norm(aj_2_r(:, end) - aj_2_r(:, r), 1)/norm(aj_2_r(:, end), 1);
        
    end
    
    %%% plotting iteraitve solution
    
    
    figure(2*j-1)
    plot(x_col1/L1, real(aj_1_r(:, end)), ':*', 'DisplayName', 'quasi true solution interpolate')
    title(['Solution on $\Gamma_{1}$, approx with ', num2str(1/C_wl(j)), 'dof per wl'])
    legend show
%     plot(x_sample1/L1, real(aj_1_r_true_it(:, end)), ':*', 'DisplayName', 'quasi true solution')
    
   
    figure(2*j)
    title(['Solution on $\Gamma_{2}$, approx with ', num2str(1/C_wl(j)), 'dof per wl'])
     plot(x_col2/L2, real(aj_2_r(:, end)), ':*', 'DisplayName', 'quasi true solution interpolate')
     legend show
end


err_normG1_it_interpl1
err_normG2_it_interpl1
%%
for j = 1:length(C_wl)-1
    
    for r = 1:R_it
        
        EOC_G1_wrt_r(j, r) = log2(err_normG1_it_interpl1(r, j)/err_normG1_it_interpl1(r, j+1));
        
        EOC_G2_wrt_r(j, r) = log2(err_normG2_it_interpl1(r, j)/err_normG2_it_interpl1(r, j+1));

    end
end

EOC_G1_wrt_r

EOC_G2_wrt_r

cond_A_G1

cond_A_G2

cond_A_S21

cond_A_S12


%% plotting the errors
x_G1_err_plots = [0:2:2*R_it-1];
figure()
hold on
for j = 1:length(C_wl)
    
    semilogy(x_G1_err_plots, err_normG1_it_interpl1(:, j), 'LineStyle', '-.', 'DisplayName', sprintf(' dof = %g', 1/C_wl(j)) )
    
    title('Approximation of $\phi_{1}^{(r)}$ on $\Gamma_{1}$ to the "true solution" with the same degrees of freedom as the approximation where r = 58')
end
set(gca, 'YScale', 'log') % But you can explicitly force it to be logarithmic
xlabel('r')
ylabel('Normalised $\ell_{1}$ error')
xlim([-0.1 x_G1_err_plots(end)+0.1])
ylim([-0.1 0.8])
legend show

x_G2_err_plots = [1:2:2*R_it];
figure()
hold on
for j = 1:length(C_wl)
    
    semilogy(x_G2_err_plots, err_normG2_it_interpl1(:, j), 'LineStyle', '-.', 'DisplayName', sprintf(' dof = %g', 1/C_wl(j)) )
    
    title('Approximation of $\phi_{2}^{(r)}$ on $\Gamma_{2}$ to the "true solution" with the same degrees of freedom as the approximation where r = 59 ')
end
set(gca, 'YScale', 'log') % But you can explicitly force it to be logarithmic
% grid on
xlabel('r')
ylabel('Normalised $\ell_{1}$ error')
xlim([-0.1 x_G2_err_plots(end)+0.1])
ylim([-0.0001 0.02])
legend show

figure()
for j = 1:length(C_wl)
    subplot(2, 1, 1)
    hold on
    semilogy(x_G1_err_plots, err_normG1_it_interpl1(:, j), 'LineStyle', '-.', 'DisplayName', sprintf(' dof = %g', 1/C_wl(j)) )

    subplot(2, 1, 2)
    hold on
    semilogy(x_G2_err_plots, err_normG2_it_interpl1(:, j), 'LineStyle', '-.', 'DisplayName', sprintf(' dof = %g', 1/C_wl(j)) )

    
end
subplot(2, 1, 1)
set(gca, 'YScale', 'log') % But you can explicitly force it to be logarithmic
xlabel('r')
ylabel('Normalised $\ell_{1}$ error')
xlim([-0.5 x_G1_err_plots(end)+0.1])
ylim([1e-11 1])
legend show
ax = gca; 
ax.FontSize = 16; 

subplot(2, 1, 2)
set(gca, 'YScale', 'log') % But you can explicitly force it to be logarithmic
xlabel('r')
ylabel('Normalised $\ell_{1}$ error')
xlim([0 x_G2_err_plots(end)+0.1])
ylim([1e-11 1e-1])
legend show
ax = gca; 
ax.FontSize = 16; 

figure(100)
xlim([-0.1 1.1])
legend show

figure(101)
xlim([-0.1 1.1])
legend show


%% Plot for internoise paper, comparison with respect to r on bndy
% In this section have a plot that looks at a fixed wavelength and plots 
% the solution with respect to the iterations, with a limited number of 
% reflections, best to also add in the "true solution".

figure()
subplot(2, 1, 1)
xlabel('$s/L_{1}$')
ylabel('$Re(\phi_{1}^{(r)}(s))$')
hold on
subplot(2, 1, 2)
xlabel('$s/L_{2}$')
ylabel('$Re(\phi_{2}^{(r)}(s))$')
hold on

% plotting first 5 solutions on the boundary
for j = 1:5
   subplot(2, 1, 1)
   plot(x_col1/L1, real(aj_1_r(:, j)), 'DisplayName', sprintf('r = %g', 2*j - 2))
   subplot(2, 1, 2)
   plot(x_col2/L2, real(aj_2_r(:, j)), 'DisplayName', sprintf('r = %g', 2*j - 1))
end
subplot(2, 1, 1)
plot(x_col1/L1, real(aj_1_r(:, end)), 'DisplayName', 'r = 58', 'LineStyle', '-.', 'Color', 'k')
xlim([-0.1 1.1])
legend show
ax = gca; 
ax.FontSize = 16; 
subplot(2, 1, 2)
plot(x_col2/L2, real(aj_2_r(:, end)), 'DisplayName', 'r = 59', 'LineStyle', '-.', 'Color', 'k')
xlim([-0.1 1.1])
legend show
ax = gca; 
ax.FontSize = 16; 
