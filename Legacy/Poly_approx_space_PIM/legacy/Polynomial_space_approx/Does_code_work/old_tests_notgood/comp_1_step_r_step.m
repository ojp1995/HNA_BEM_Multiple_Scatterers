%%% Test 2, here we will be comparing the iterative method to the method
%%% where we solve it in one
clear all

% loading the true solution:
load('SC_alin1_phi1_phi2_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_m1')

addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx')

theta = 0; %OC angle between the downwards verticle and incident wave, anti clockwise
d_test = [sin(theta), - cos(theta)];  % simon code, direction of incident wave

if d_test == d
    disp('direction of incident wave matches')
else
    error('d loaded does not match theta for thos case')
end
clear C_wl

C_wl = [1/5, 1/10, 1/20, 1/40, 1/80, 1/160]; % h <= C_wl \lambda, 1/C_wl intervals per wl.

% introducing the problem parmeters:
G1 = [r1start, r1end ];
G2 = [r2start, r2end ];

if L1==sqrt( (G1(3) - G1(1))^2 + (G1(4) - G1(2))^2 )
else
    error('G1 does not match data given')
end

if L2==sqrt( (G2(3) - G2(1))^2 + (G2(4) - G2(2))^2 )
else
    error('G2 does not match data given')
end


phi1_true = phi1;
phi2_true = phi2;

% k = 5;  % wavenumber = 1, as if we have nondimensionalised


% constants needed for the smoothing function
C1 = 1;
C2 = pi;



R = 3; % this is the number of re-reflections we want


N_sample = length(phi1) - 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% All on 1 solve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading true solution!
% load()
% set up for simons code:

% % % % r1start = [G1(1),G1(2)]; r1end = [G1(3),G1(4)]; % Endpoints of screen 1
% % % % r2start = [G2(1),G2(2)]; r2end = [G2(3),G2(4)]; %Endpoints of screen 2
% % % % d = [sin(theta),-cos(theta)]; % the direction of the incident wave
% % % % N1 = 8000; N2 = 8000; % the number of boundary elements on screens 1 and 2
% % % % h1 = norm(r1end-r1start)/N1 % the length of each element on screen 1
% % % % h2 = norm(r2end-r2start)/N2 % the length of each element on screen 2
% % % % Next calculate the x and y coordinates of the element midpoints on screen
% % % % 1 and then on screen 2
% % % % x1 = r1start(1)+((1:N1)-0.5)*(r1end(1)-r1start(1))/N1;
% % % % y1 = r1start(2)+((1:N1)-0.5)*(r1end(2)-r1start(2))/N1;
% % % % x2 = r2start(1)+((1:N2)-0.5)*(r2end(1)-r2start(1))/N2;
% % % % y2 = r2start(2)+((1:N2)-0.5)*(r2end(2)-r2start(2))/N2;
% % % % h = [h1*ones(size(x1)),h2*ones(size(x2))]; % the lengths of the elements
% % % % x = [x1,x2]; % the x-coords of the element midpoints
% % % % y = [y1,y2]; % and their y-coords
% % % % phi = bem2(x,y,h,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
% % % % phi1 = phi(1:N1); phi2 = phi(N1+1:end); % Split phi into the solution on screens 1 and 2
% % % % isolating solutions:

 



% OP code, in theory this should converge well at every point
% A = [A1 S12; S21 A2];
% u_inc = zeros(length(N1+N2), 1);
% u_inc(1:N1, 1) = PW_incident(k, theta, G1, x_col1);
% u_inc(N1+1:N1+N2, 1) = PW_incident(k, theta, G2, x_col2);
% 
% aj_12 = A\u_inc;
% 
% aj_1 = aj_12(1:N1, 1);
% aj_2 = aj_12(N1+1:N1+N2, 1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:length(C_wl)
    r = 1;
    % number of intervals on wach screen:
     N1 = ceil(k*L1./(C_wl(j)*2*pi)) % number of itervals on G1
     N2 = ceil(k*L2./(C_wl(j)*2*pi)) % number of intervals on G2



    Q = 1;  % quadrature discretisation

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
    % phi1_r = zeros(N_sample, R);
    % phi2_r = zeros(N_sample, R);
    % r = 1 case:
    u_i_1_r(:, 1) = PW_incident(k, theta, G1, x_col1);
    aj_1_r(:, 1) = A1\u_i_1_r(:, 1);
    [phi1_r(:, 1), x_sample1] = coeff_2soln_midpoint(aj_1_r(:, 1), L1, N_sample, N1);


    u_i_2_r(:, 1) =  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_r(:, 1), L1, N1-1, N1);
    aj_2_r(:, 1) = A2\u_i_2_r(:, 1);
    [phi2_r(:, 1), x_sample2] = coeff_2soln_midpoint(aj_2_r(:, 1), L2, N_sample, N2);

    phi1 = interp1(x_sample1, phi1_true, x_col1);
    phi2 = interp1(x_sample2, phi2_true, x_col2);
    
    
    % computing errors:
%     err_G1(j, 1) = norm(abs( phi1 - aj_1_r(:, 1) ), 2);
%     err_G2(j, 1) = norm(abs( phi2 - phi2_r(:, 1) ), 2);
% 
%     err_norm_G1(j, 1) = norm(abs( phi1(1:end-1) - aj_1_r(1:end-1, 1) ), 2)/norm(phi1(1:end-1), 2);
%     err_norm_G2(j, 1) = norm(abs( phi2(1:end-1) - phi2_r(1:end-1, 1) ), 2)/norm(phi2(1:end-1), 2);
%     
    
    err_norm_G12_interpB2S(j, 2*r-1) = norm( phi1 - aj_1_r(:, 1).' , 1)/norm(phi1, 1);
    err_norm_G12_interpB2S(j, 2*r) = norm(phi2 - aj_2_r(:, 1).' , 1)/norm(phi2, 1);
    
    % does not work, intper leads to a load of NaNs
%     phi1_interp = interp1(x_col1, aj_1_r(:, r), x_sample1);
%     phi2_interp = interp1(x_col2 , aj_2_r(:, r), x_sample2 );
%     
%     err_norm_G12_interpS2B(j, 2*r-1) = norm(phi1_true -  phi1_interp.', 1)/norm(phi1_true, 1);
%     err_norm_G12_interpS2B(j, 2*r) = norm( phi2_true - phi2_interp.', 1)/norm(phi2_true, 1);

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

    for r = 2:R
    %     compute the even solutions on \Gamma_{1}
    %     first compute the incident:
        u_i_1_r(:, r) =  PW_incident(k, theta, G1, x_col1).' - S12*coeff_2soln_midpoint(aj_2_r(:, r-1), L2, N2-1, N2);
    %     conpute the coefficients
        aj_1_r(:, r) = A1\u_i_1_r(:, r);
    %     compute the solution
         [phi1_r(:, r), ~] = coeff_2soln_midpoint(aj_1_r(:, r), L1, N_sample, N1);
%         figure(5)
%         plot(x_col1/L1, aj_1_r(:, r), 'DisplayName', sprintf('r = %g', 2*r - 2))

    %     compute the odd solutions on \Gamma_{2}
    %         first compute the incident:
        u_i_2_r(:, r) =  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_r(:, r), L1, N1-1, N1);
    %     conpute the coefficients
        aj_2_r(:, r) = A2\u_i_2_r(:, r);
    %     compute the solution
         [phi2_r(:, r), ~] = coeff_2soln_midpoint(aj_2_r(:, r), L2, N_sample, N2);

        % computing errors:
%         err_G1(j, r) = norm(abs( phi1 - phi1_r(:, r) ), 2);
%         err_G2(j, r) = norm(abs( phi2 - phi2_r(:, r) ), 2);
% 
%         err_norm_G1(j, r) = norm(abs( phi1 - phi1_r(:, r) ), 2)/norm(phi1, 2);
%         err_norm_G2(j, r) = norm(abs( phi2 - phi2_r(:, r) ), 2)/norm(phi2, 2);
        
        disp(r)
        
        %%%% Could this help? Would seem that last value of OC phi approx
        %%%% is 0 so this removes it.
%         err_G1(j, r) = norm(abs( phi1(1:end-1) - phi1_r(1:end-1, r) ), 2);
%         err_G2(j, r) = norm(abs( phi2(1:end-1) - phi2_r(1:end-1, r) ), 2);
% 
%         err_norm_G1(j, r) = norm(abs( phi1(1:end-1) - phi1_r(1:end-1, r) ), 2)/norm(phi1(1:end-1), 2);
%         err_norm_G2(j, r) = norm(abs( phi2(1:end-1) - phi2_r(1:end-1, r) ), 2)/norm(phi2(1:end-1), 2);
%         
        
        phi1 = interp1(x_sample1, phi1_true, x_col1);
        phi2 = interp1(x_sample2, phi2_true, x_col2);
        
        err_norm_G12_interpB2S(j, 2*r-1) = norm(phi1 - aj_1_r(:, r).' , 1)/norm(phi1, 1);
        err_norm_G12_interpB2S(j, 2*r) = norm( phi2 - aj_2_r(:, r).' , 1)/norm(phi2, 1);
        
%         phi1_interp = interp1(x_col1, aj_1_r(:, r), x_sample1);
%         phi2_interp = interp1(x_col2 , aj_2_r(:, r), x_sample2 );

%         err_norm_G12_interpS2B(j, 2*r-1) = norm(phi1_true -  phi1_interp, 1)/norm(phi1, 1);
%         err_norm_G12_interpS2B(j, 2*r) = norm( phi2_true - phi2_interp , 1)/norm(phi2, 1);
%          figure(6)
%         plot(x_col2/L2, aj_2_r(:, r), 'DisplayName', sprintf('r = %g', 2*r-1))
% 
%         figure(7)
%         subplot(2, 1, 1)
%         plot(x_col1/L1, aj_1_r(:, r), 'DisplayName', sprintf('r = %g', 2*r - 2))
% 
%         subplot(2, 1, 2)
%         plot(x_col2/L2, aj_2_r(:, r), 'DisplayName', sprintf('r = %g', 2*r-1))

    end

%     figure(5)
%     plot(x_col1/L1, aj_1, 'DisplayName', 'Single solve')
%     legend show
%     title('Solution on \Gamma_{1} for increasing number of iterations')
%     xlabel('x')
%     ylabel('\phi_{1}^{r}')
%     xlim([-0.1 1.1])
%     figure(6)
%     plot(x_col2/L2, aj_2, 'DisplayName', 'Single solve')
%     title('Solution on \Gamma_{2} for increasing number of iterations')
%     xlabel('x')
%     ylabel('\phi_{2}^{r}')
%     xlim([-0.1 1.1])
%     legend show
% 
%     figure(7)
%     subplot(2, 1, 1)
%     title('Solution on \Gamma_{1} for increasing number of iterations')
%     xlabel('$s/L_{1}$', 'FontSize', 17, 'interpreter', 'latex')
%     ylabel('Re $\left(\phi_{1}^{r} \right)$', 'FontSize', 17, 'interpreter', 'latex')
%     xlim([-0.1 1.1])
%     legend show
%     set(legend,'fontsize',17);
%     subplot(2, 1, 2)
%     title('Solution on \Gamma_{2} for increasing number of iterations')
%     xlabel('$s/L_{2}$', 'FontSize', 17, 'interpreter', 'latex')
%     ylabel('Re $ \left(\phi_{2}^{r} \right) $', 'FontSize', 17, 'interpreter', 'latex')
%     xlim([-0.1 1.1])
%     legend show
%     set(legend,'fontsize',17);
end
% err_G1
% err_G2
% 
% err_norm_G1
% err_norm_G2

err_norm_G12_interpB2S

% now comuting EOC
for r = 1:2*R
    for j = 1:(length(C_wl) - 1)
        
        EOC_interpB2S(j, r) = log2(err_norm_G12_interpB2S(j+1, r))/log2(err_norm_G12_interpB2S(j, r));
        
    end
    
end
