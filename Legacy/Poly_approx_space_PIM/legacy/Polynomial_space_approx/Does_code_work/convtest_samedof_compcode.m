% in this script we will be testing the convergence between my code and
% simons code for each iteration. The dof will be the same for each, the
% reason being as this is using the PIM then it is only valid (might not be
% the right word) at the midpoints.

clear all

%loading paths
addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/') % code to be tested
addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/Simon_code/') % simons code, taken to be as working code

% introducing the problem parameters:

k = 5; % wavenumber
theta = 0; %OC angle between the downwards verticle and incident wave, anti clockwise
d = [sin(theta), - cos(theta)];  % simon code, direction of incident wave

% range of dof we want to loop over
C_wl = [1/5, 1/10, 1/20, 1/40, 1/80, 1/160];  % length of each interval C_wl a wavelength long
% C_wl = [1/5, 1/10];
%%%%% Orientation of the screen
%SC
r1start = [-2*pi, 2*pi]; r1end = [0, 0]; % \Gamma_{1}
r2start = [7*pi, 3*pi]; r2end = [4*pi, 0]; % \Gamma_{2}


% OC
G1 = [r1start, r1end ];
G2 = [r2start, r2end ];

L1 = sqrt( (G1(3) - G1(1))^2 + (G1(4) - G1(2))^2 );  % length of G1
L2 = sqrt( (G2(3) - G2(1))^2 + (G2(4) - G2(2))^2 );  % length of G2

C1 = 1;
C2 = pi;

% now we can begin the looping over the range of dofs per wl
for j = 1:length(C_wl)
    
    % discretisation:
    N1 = ceil(k*L1./(C_wl(j)*2*pi)) % number of itervals on G1
    N2 = ceil(k*L2./(C_wl(j)*2*pi)) % number of intervals on G2
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%  Simon code   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h1 = norm(r1end-r1start)/N1; % the length of each element on screen 1
    h2 = norm(r2end-r2start)/N2; % the length of each element on screen 2
%     Next calculate the x and y coordinates of the element midpoints on screen
%     1 and then on screen 2
    
    x1 = r1start(1)+((1:N1)-0.5)*(r1end(1)-r1start(1))/N1;
    y1 = r1start(2)+((1:N1)-0.5)*(r1end(2)-r1start(2))/N1;
    x2 = r2start(1)+((1:N2)-0.5)*(r2end(1)-r2start(1))/N2;
    y2 = r2start(2)+((1:N2)-0.5)*(r2end(2)-r2start(2))/N2;
    h1vector = h1*ones(size(x1)); % Lengths of the elements on screen 1
    h2vector = h2*ones(size(x2)); % and on screen 2
    h = [h1vector,h2vector]; % the lengths of the elements
    x = [x1,x2]; % the x-coords of the element midpoints
    y = [y1,y2]; % and their y-coords
    
%     step 0:
    phitrue1_0 = bem2(x1,y1,h1vector,k,d);
    
%     step 1:
    phitrue2_1 = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phitrue1_0,k,d); % approximations to the values of phi at the element midppoints, phi(j) the value at (x(j),y(j))
    
%     step 2:
    phitrue1_2 = bem2BS(x2,y2,x1,y1,h2vector,h1vector,phitrue2_1,k,d);
    
%     step 3:
    phitrue2_3 = bem2BS(x1,y1,x2,y2,h1vector,h2vector,phitrue1_2,k,d);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%  Oliver code   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Q = 1;
    [A1, x_col1, A1_sing, A1_smooth] = screen_mat_poly_approx_PIM(N1, 0, L1, k, C1, C2, Q);
    [A2, x_col2, A2_sing, A2_smooth] = screen_mat_poly_approx_PIM(N2, 0, L2, k, C1, C2, Q);
    
    [S21, col_points_1] = S21_op(G1, G2, k, N1, N2);
    [S12, col_points_2] = S21_op(G2, G1, k, N2, N1);
    
    
    % step 0
    u_i_1_0 = PW_incident(k, theta, G1, x_col1);
    aj_1_0 = A1\u_i_1_0.';
%     [phi_1_0, x1_plot] = coeff_2soln_midpoint(aj_1_0, L1, N_sample, N1);
    
    % step 1
    u_i_2_1 =  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_0, L1, N1-1, N1);
    aj_2_1 = A2\u_i_2_1;
%     [phi_2_1, x2_plot] = coeff_2soln_midpoint(aj_2_1, L2, N_sample, N2);

    % step 2
    u_i_1_2 =  PW_incident(k, theta, G1, x_col1).' - S12*coeff_2soln_midpoint(aj_2_1, L2, N2-1, N2);
    aj_1_2 = A1\u_i_1_2;
%     phi_1_2 = coeff_2soln_midpoint(aj_1_2, L1, N_sample, N1);
    
    % step 3
    u_i_2_3=  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_2, L1, N1-1, N1);
    aj_2_3= A2\u_i_2_3;
%     phi_2_3 = coeff_2soln_midpoint(aj_2_3, L2, N_sample, N2);
    

    % Now computing errors between the two
    err_norml1(j, 1) = norm(abs(phitrue1_0 - aj_1_0), 1)/norm(phitrue1_0, 1);
    err_norml1(j, 2) = norm(abs(phitrue2_1 - aj_2_1), 1)/norm(phitrue2_1, 1);
    err_norml1(j, 3) = norm(abs(phitrue1_2 - aj_1_2), 1)/norm(phitrue1_2, 1);
    err_norml1(j, 4) = norm(abs(phitrue2_3 - aj_2_3), 1)/norm(phitrue2_3, 1);
    
        % Now computing errors between the two
    err_norml2(j, 1) = norm(abs(phitrue1_0 - aj_1_0), 2)/norm(phitrue1_0, 2);
    err_norml2(j, 2) = norm(abs(phitrue2_1 - aj_2_1), 2)/norm(phitrue2_1, 2);
    err_norml2(j, 3) = norm(abs(phitrue1_2 - aj_1_2), 2)/norm(phitrue1_2, 2);
    err_norml2(j, 4) = norm(abs(phitrue2_3 - aj_2_3), 2)/norm(phitrue2_3, 2);
    
           % Now computing errors between the two
    err_normlinf(j, 1) = norm(abs(phitrue1_0 - aj_1_0), inf)/norm(phitrue1_0, inf);
    err_normlinf(j, 2) = norm(abs(phitrue2_1 - aj_2_1), inf)/norm(phitrue2_1, inf);
    err_normlinf(j, 3) = norm(abs(phitrue1_2 - aj_1_2), inf)/norm(phitrue1_2, inf);
    err_normlinf(j, 4) = norm(abs(phitrue2_3 - aj_2_3), inf)/norm(phitrue2_3, inf);
   
%     plotting
    figure(1)
    plot(x_col1/L1, real(aj_1_0), 'DisplayName', sprintf('OC dof= %g', 1/C_wl(j)))
    hold on
    plot(x_col1/L1, real(phitrue1_0), 'DisplayName', sprintf('SC dof= %g', 1/C_wl(j)))
    figure(2)
    plot(x_col2/L2, real(aj_2_1), 'DisplayName', sprintf('OC dof= %g', 1/C_wl(j)))
    hold on
    plot(x_col2/L2, real(phitrue2_1), 'DisplayName', sprintf('SC dof= %g', 1/C_wl(j)))
    
    figure(3)
    plot(x_col1/L1, real(aj_1_2), 'DisplayName', sprintf('OC dof= %g', 1/C_wl(j)))
    hold on
    plot(x_col1/L1, real(phitrue1_2), 'DisplayName', sprintf('SC dof= %g', 1/C_wl(j)))
    figure(4)
    plot(x_col2/L2, real(aj_2_3), 'DisplayName', sprintf('OC dof= %g', 1/C_wl(j)))
    hold on
    plot(x_col2/L2, real(phitrue2_3), 'DisplayName', sprintf('SC dof= %g', 1/C_wl(j)))
    
end
figure(1)
title('Approximations for $\phi_{1}^{(0)}$')
legend show

figure(2)
title('Approximations for $\phi_{2}^{(1)}$')
legend show

figure(3)
title('Approximations for $\phi_{1}^{(2)}$')
legend show

figure(4)
title('Approximations for $\phi_{2}^{(3)}$')
legend show
R = 4;
% computing EOC
for j = 1:(length(C_wl)-1)
    for r = 1:R
        EOC_l1(j, r) = log2(err_norml1(j+1, r)/err_norml1(j, r));
        EOC_l2(j, r) = log2(err_norml2(j+1, r)/err_norml2(j, r));
        EOC_linf(j, r) = log2(err_normlinf(j+1, r)/err_normlinf(j, r));
    end
end

err_norml1
EOC_l1
err_norml2
EOC_l2
err_normlinf
EOC_linf

