% this script is to compute the scattering problem in a polynomial
% approximation space, approxmating the integrals using a product
% integration method and the mid point rule.

clear all

% figure(1)
% figure(2)
% figure(3)
% 
% close figure 1
% close figure 2
% close figure 3

tic
% Introducing the screens
% G1 = [-0.1, 2*pi, 0, 0];
G1 = [-2*pi, 2*pi, 0, 0];
% G1 = [-7*pi, 8*pi, 1*pi, 0];
% G1 = [-8*pi, 8*pi, 0, 0];
% G1 = [0, 0, 2*pi, 0];  % coordinates for \Gamma_{1
% G1 = [-8*pi, 8*pi, 0, 0];
% G1 = [0, 8*pi, 1*pi, 0];
% G1 = [0, 0, 2*pi, 0];
L1 = sqrt( (G1(3) - G1(1))^2 + (G1(4) - G1(2))^2 );  % length of G1

% N1 = 50;  % number of intervals on G1
% h1 = L1/N1;  % step size on G1

d = 4*pi; % distance between the screens
% G2 = [2*pi, 0, 2*pi+0.5, 3*pi];
G2 = [ 2*pi, 0, 5*pi, 3*pi];
% G2 = [d, 0, d+12*pi, 12*pi];
% G2 = [8*pi, 8*pi, 16*pi, 0];
% G2 = [d, 0, d + 8*pi, 0];  % coordinates for \Gamma_{2}
% G2 = [2*pi, 10*pi, 11*pi, 10/*pi];
% G2 = [4*pi, 0, 5*pi, 8*pi];
% G2 = [pi/4, 2*pi, 3*pi, 2*pi];
L2 = sqrt( (G2(3) - G2(1))^2 + (G2(4) - G2(2))^2 );  % length of G2

C_wl= 1/20

k = 10;  % wavenumber


N1 = ceil(k*L1./(C_wl*2*pi)) % number of itervals on G1
N2 = ceil(k*L2./(C_wl*2*pi)) % number of intervals on G2
    

% N2 = N1;  % number of intervals on G2
% % h2 = L2/N2;  % step size on G2


theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

% I think the set up has now been introduced so we can start to compute the
% matrices for \Gamma_{1} and \Gamma_{2}

Q = 1;  % quadrature discretisation

[A1, x_col1 A1_sing, A1_smooth] = screen_mat_poly_approx_PIM(N1, 0, L1, k, C1, C2, Q);
[A2, x_col2, A2_sing, A2_smooth] = screen_mat_poly_approx_PIM(N2, 0, L2, k, C1, C2, Q);


% figure()
% surf(A1_sing)
% shading interp
% 
% figure()
% surf(real(A1_smooth))
% shading interp
% computing the S21 and S12 operators
[S21, col_points_1] = S21_op(G1, G2, k, N1, N2);
[S12, col_points_2] = S21_op(G2, G1, k, N2, N1);
% % Now lets compute a right hand side
% u_i_1_0 =  PW_incident(k, theta, G1, x_col1);
% % solve for \Gamma_{1}
% aj1_0 = A1\u_i_1_0.';
% 
% % first re-reflection
% u_i_2_1 =  PW_incident(k, theta, G2, x_col1).' - S21*coeff_2soln_midpoint(aj1_0, L1, N2-1, N1);
% aj2_1 = A2\u_i_2_1;
% 
% % second reflection
% u_i_1_2 = PW_incident(k, theta, G1, x_col1).' - S12*coeff_2soln_midpoint(aj2_1, L2, N1-1, N2);
% aj1_2 = A1\u_i_1_2;
% 
% % third reflection
% u_i_2_3 =  PW_incident(k, theta, G2, x_col1).' - S21*coeff_2soln_midpoint(aj1_2, L1, N2-1, N1);
% aj2_3 = A2\u_i_2_3;
% 
% % fourth reflection
% u_i_1_4 = PW_incident(k, theta, G1, x_col1).' - S12*coeff_2soln_midpoint(aj2_3, L2, N1-1, N2);
% aj1_4 = A1\u_i_1_4;
% 
% % fifth reflection
% u_i_2_5 =  PW_incident(k, theta, G2, x_col1).' - S21*coeff_2soln_midpoint(aj1_4, L1, N2-1, N1);
% aj2_5 = A2\u_i_2_5;




% compute the solutions from coeff2soln
N_sample = 1000;  % number of points we are sampling the solution at

% [phi1_0, x_sample0] = coeff_2soln_midpoint(aj1_0, L1, N_sample, N1);
% [phi2_1, x_sample1] = coeff_2soln_midpoint(aj2_1, L2, N_sample, N2);
% [phi1_2, x_sample2] = coeff_2soln_midpoint(aj1_2, L1, N_sample, N1);
% [phi2_3, x_sample3] = coeff_2soln_midpoint(aj2_3, L2, N_sample, N2);
% [phi1_4, x_sample4] = coeff_2soln_midpoint(aj1_4, L1, N_sample, N1);
% [phi2_5, x_sample5] = coeff_2soln_midpoint(aj2_5, L2, N_sample, N2);
% [phi1_6, x_sample6] = coeff_2soln_midpoint(aj1_2, L1, N_sample, N1);
% toc

% automating the re-reflections
R = 4; % this is the number of re-reflections we want
% initialising
u_i_1_r = zeros(length(x_col1), R);
aj_1_r = zeros(length(x_col1), R);
u_i_2_r = zeros(length(x_col2), R);
aj_2_r = zeros(length(x_col2), R);
phi1_r = zeros(N_sample, R);
phi2_r = zeros(N_sample, R);
% r = 1 case:
u_i_1_r(:, 1) = PW_incident(k, theta, G1, x_col1);
aj_1_r(:, 1) = A1\u_i_1_r(:, 1);
[phi1_r(:, 1), x_sample1] = coeff_2soln_midpoint(aj_1_r(:, 1), L1, N_sample-1, N1);


u_i_2_r(:, 1) =  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_r(:, 1), L1, N1-1, N1);
aj_2_r(:, 1) = A2\u_i_2_r(:, 1);
[phi2_r(:, 1), x_sample2] = coeff_2soln_midpoint(aj_2_r(:, 1), L2, N_sample-1, N2);

%%
figure(5);
plot(x_col1/L1, real(aj_1_r(:, 1)), 'DisplayName', 'r = 0')
hold on


figure(6);
plot(x_col2/L2, real(aj_2_r(:, 1)), 'DisplayName', 'r = 1')
hold on

figure(7)
subplot(2, 1, 1)
plot(x_col1/L1, aj_1_r(:, 1), 'DisplayName', sprintf('r = %g', 0))
hold on
subplot(2, 1, 2)
plot(x_col2/L2, aj_2_r(:, 1), 'DisplayName', sprintf('r = %g', 1))
hold on
% 
for r = 2:R
%     compute the even solutions on \Gamma_{1}
%     first compute the incident:
    u_i_1_r(:, r) =  PW_incident(k, theta, G1, x_col1).' - S12*coeff_2soln_midpoint(aj_2_r(:, r-1), L2, N2-1, N2);
%     conpute the coefficients
    aj_1_r(:, r) = A1\u_i_1_r(:, r);
%     compute the solution
    [phi1_r(:, r), ~] = coeff_2soln_midpoint(aj_1_r(:, r), L1, N_sample-1, N1);
    figure(5)
    plot(x_col1/L1, aj_1_r(:, r), 'DisplayName', sprintf('r = %g', 2*r - 2))
    
%     compute the odd solutions on \Gamma_{2}
%         first compute the incident:
    u_i_2_r(:, r) =  PW_incident(k, theta, G2, x_col2).' - S21*coeff_2soln_midpoint(aj_1_r(:, r), L1, N1-1, N1);
%     conpute the coefficients
    aj_2_r(:, r) = A2\u_i_2_r(:, r);
%     compute the solution
    [phi2_r(:, r), ~] = coeff_2soln_midpoint(aj_2_r(:, r), L2, N_sample-1, N2);
    
     figure(6)
    plot(x_col2/L2, aj_2_r(:, r), 'DisplayName', sprintf('r = %g', 2*r-1))
    
    figure(7)
    subplot(2, 1, 1)
    plot(x_col1/L1, aj_1_r(:, r), 'DisplayName', sprintf('r = %g', 2*r - 2))
    
    subplot(2, 1, 2)
    plot(x_col2/L2, aj_2_r(:, r), 'DisplayName', sprintf('r = %g', 2*r-1))
    
end


figure(5)
legend show
title('Solution on \Gamma_{1} for increasing number of iterations')
xlabel('x')
ylabel('\phi_{1}^{r}')
xlim([-0.1 1.1])
figure(6)
title('Solution on \Gamma_{2} for increasing number of iterations')
xlabel('x')
ylabel('\phi_{2}^{r}')
xlim([-0.1 1.1])
legend show

figure(7)
subplot(2, 1, 1)
% title('Solution on \Gamma_{1} for increasing number of iterations')
xlabel('$s/L_{1}$', 'FontSize', 17, 'interpreter', 'latex')
ylabel('Re $\left(\phi_{1}^{r}(s) \right)$', 'FontSize', 17, 'interpreter', 'latex')
xlim([-0.1 1.1])
legend show
set(legend,'fontsize',17);
subplot(2, 1, 2)
% title('Solution on \Gamma_{2} for increasing number of iterations')
xlabel('$s/L_{2}$', 'FontSize', 17, 'interpreter', 'latex')
ylabel('Re $ \left(\phi_{2}^{r}(s) \right) $', 'FontSize', 17, 'interpreter', 'latex')
xlim([-0.1 1.1])
legend show
set(legend,'fontsize',17);
% figure(); plot(x_sample0, real(phi1_0), x_sample2, real(phi1_2), x_sample4, real(phi1_4)); legend('\phi_{1}^{(0)}', '\phi_{1}^{(2)}', '\phi_{1}^{(3)}');legend show; title('Solutions on \Gamma_{1}')
% figure(); plot(x_sample1, real(phi2_1), x_sample3, real(phi2_3), x_sample5, real(phi2_5)); legend('\phi_{2}^{(1)}', '\phi_{2}^{(3)}', '\phi_{2}^{(5)}');legend show; title('Solutions on \Gamma_{2}')

%%
keyboard
% compute the solution in the domain
% Makr a square grid around the screens
X_coordinates = [G1(1), G1(3), G2(1), G2(3)];
X_min = min(X_coordinates);
X_max = max(X_coordinates);
X_diff = abs(X_max - X_min);

Y_coordinates = [G1(2), G1(4), G2(2), G2(4)];
Y_min = min(Y_coordinates);
Y_max = max(Y_coordinates);
Y_diff = abs(Y_max - Y_min);

% maybe also change this part
N_diff = max(abs(Y_diff - X_diff)) + 5*pi;  % largest distance either x direction or y + 4*pi (2*pi in either direction)

N_d = 500;

min_leftorbottom = min(X_min, Y_min);
max_toporright = max(X_max, Y_max);
h_d = N_diff/N_d;
% need to change this bit so that it is closer zoomed in on the picture
X = [min_leftorbottom - 5: h_d: max_toporright + 5];
Y = [min_leftorbottom - 5: h_d: max_toporright + 5];
% first compute the incident field

U_i = incident_domain(k, theta, X, Y);

figure(); pcolor(X, Y, real(U_i)); shading interp



toc

% % % tic
% % % computing solution in 1 step:
%%% OP code, in theory this should converge well at every point
% % % A = [A1 S12; S21 A2];
% % % u_inc = zeros(length(N1+N2), 1);
% % % u_inc(1:N1, 1) = PW_incident(k, theta, G1, x_col1);
% % % u_inc(N1+1:N1+N2, 1) = PW_incident(k, theta, G2, x_col2);
% % % 
% % % aj_12 = A\u_inc;
% % % 
% % % aj_1 = aj_12(1:N1, 1);
% % % aj_2 = aj_12(N1+1:N1+N2, 1);
% % % 
% % % toc

% compute the scattered field
Q_d = 110;



% NEWOPC_phi12_itall_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_n1 = aj_1_r;
% 
% save('NEWOPC_phi12_itall_g1_n2pi_2pi_0_0_G2_4pi_0_7pi_3pi_Cwl640_d0_n1')
% 

tic
[us1_0, x1, y1] = scattered_domain( G1, aj_1_r(:, 1), k, X, Y, N1);
figure(); pcolor(X, Y, real(us1_0)); shading interp; colorbar; title('scattered field no-rereflections ')
toc

tic
% [us2_1, x2, y2] = scattered_domain( G2, aj_2_r(:, 1), k, X, Y, N2);
[us2_1, x1, y1, x2, y2] = scattered_domain2( G1, G2, aj_1_r(:, 1), aj_2_r(:, 1), k, X, Y, N1, N2);
figure(); pcolor(X, Y, real(us2_1)); shading interp; colorbar; title('scattered field 1 rereflections ')
toc

tic
us1_2 = scattered_domain( G1, aj_1_r(:, 2), k, X, Y, N1);
[us1_2, x1, y1, x2, y2] = scattered_domain2( G1, G2, aj_1_r(:, 2), aj_2_r(:, 1), k, X, Y, N1, N2);
figure(); pcolor(X, Y, real(us1_2)); shading interp; colorbar; title('scattered field 2 rereflections ')
toc

tic
us2_3 = scattered_domain2( G1, G2, aj_1_r(:, 2), aj_2_r(:, 2), k, X, Y, N1, N2);
figure(); pcolor(X, Y, real(us2_3)); shading interp; colorbar; title('scattered field 3 rereflections ')
toc

tic
us1_4 = scattered_domain2( G1, G2, aj_1_r(:, 3), aj_2_r(:, 2), k, X, Y, N1, N2);
% figure(); pcolor(X, Y, real(us1_4)); shading interp; colorbar; title('scattered field 4 rereflections ')
toc

tic
us2_5 = scattered_domain2( G1, G2, aj_1_r(:, 3), aj_2_r(:, 3), k, X, Y, N1, N2);
% figure(); pcolor(X, Y, real(us2_5)); shading interp; colorbar; title('scattered field 5 rereflections ')
toc

UT_1_0 = U_i - us1_0;

UT_2_1 = U_i - us2_1;
UT_1_2 = U_i - us1_2;
UT_2_3 = U_i - us2_3;
UT_1_4 = U_i - us1_4;
UT_2_5 = U_i - us2_5;

% UT_2_1 = UT_1_0 - us2_1;
% UT_1_2 = UT_2_1 - us1_2;
% UT_2_3 = UT_1_2 - us2_3;
% UT_1_4 = UT_2_3 - us1_4;
% UT_2_5 = UT_1_4 - us2_5;

figure(); pcolor(X, Y, real(UT_1_0)); shading interp; colormap(jet); colorbar; title('total field no-rereflections ')
hold on
Gamma_1 = plot(x1, y1);
Gamma_1.LineWidth = 4;
Gamma_1.Color = [0 0 0];

% UT_2_1 = U_i - us2_1;

figure(); pcolor(X, Y, real(UT_2_1)); shading interp; colormap(jet); colorbar; title('total field 1 re-reflection ')
hold on
Gamma_1 = plot(x1, y1);
Gamma_1.LineWidth = 4;
Gamma_1.Color = [0 0 0];
Gamma_2 = plot(x2, y2);
Gamma_2.LineWidth = 4;
Gamma_2.Color = [0 0 0];

figure(); pcolor(X, Y, real(UT_1_2)); shading interp; colormap(jet); colorbar; title('total field 2 rereflections ')
hold on
Gamma_1 = plot(x1, y1);
Gamma_1.LineWidth = 4;
Gamma_1.Color = [0 0 0];
Gamma_2 = plot(x2, y2);
Gamma_2.LineWidth = 4;
Gamma_2.Color = [0 0 0];

figure(); pcolor(X, Y, real(UT_2_3)); shading interp; colormap(jet); colorbar; title('total field 3 re-reflection ')
hold on
Gamma_1 = plot(x1, y1);
Gamma_1.LineWidth = 4;
Gamma_1.Color = [0 0 0];
Gamma_2 = plot(x2, y2);
Gamma_2.LineWidth = 4;
Gamma_2.Color = [0 0 0];

figure(); pcolor(X, Y, real(UT_1_4)); shading interp; colormap(jet); colorbar; title('total field 4 rereflections ')
hold on
Gamma_1 = plot(x1, y1);
Gamma_1.LineWidth = 4;
Gamma_1.Color = [0 0 0];
Gamma_2 = plot(x2, y2);
Gamma_2.LineWidth = 4;
Gamma_2.Color = [0 0 0];

figure(); pcolor(X, Y, real(UT_2_5)); shading interp; colormap(jet); colorbar; title('total field 5 re-reflection ')
hold on
Gamma_1 = plot(x1, y1);
Gamma_1.LineWidth = 4;
Gamma_1.Color = [0 0 0];
Gamma_2 = plot(x2, y2);
Gamma_2.LineWidth = 4;
Gamma_2.Color = [0 0 0];

% figure for iterative presentation
% iterative solution for \Gamma_{1}
% figure()
% plot(x_col1/L1, real(aj_1_r(:, 1)), 'DisplayName', 'r = 0')
% title('Solution on \Gamma_{1} for increasing number of iterations')
% xlabel('x')
% ylabel('\phi_{1}^{r}')
% xlim([-0.1 1.1])
% legend show
% hold on
% plot(x_col1/L1, real(aj_1_r(:, 2)), 'DisplayName', 'r = 2')
% plot(x_col1/L1, real(aj_1_r(:, 3)), 'DisplayName', 'r = 4')
% 
% figure()
% plot(x_col2/L2, real(aj_2_r(:, 1)), 'DisplayName', 'r = 1')
% title('Solution on \Gamma_{2} for increasing number of iterations')
% xlabel('x')
% ylabel('\phi_{2}^{r}')
% xlim([-0.1 1.1])
% legend show
% hold on
% plot(x_col2/L2, real(aj_2_r(:, 2)), 'DisplayName', 'r = 3')
% plot(x_col2/L2, real(aj_2_r(:, 3)), 'DisplayName', 'r = 5')


% figure()
% for j = 1:R
%     subplot(2, 1, 1)
%     plot(x_col1/L1, aj_1_r(:, r), 'DisplayName', sprintf('r = %g', 2*r - 2))
% 
%     subplot(2, 1, 2)
%     plot(x_col2/L2, aj_2_r(:, r), 'DisplayName', sprintf('r = %g', 2*r-1))
% end

% incident on \Gamma_{2} first
% aj_2_r0 = A2\PW_incident(k, theta, G2, x_col2).';
% 
% aj_1_r1 = A1\(PW_incident(k, theta, G2, x_col2).' - S12*coeff_2soln_midpoint(aj_2_r0, L2, N2-1, N2));
% 
% tic
% us2_0 = scattered_domain( G2, aj_2_r0, k, X, Y, N2);
% figure(); pcolor(X, Y, real(us2_0)); shading interp; colorbar; title('scattered field no-rereflections, \Gamma_{2} ')
% toc
% 
% tic
% us1_1 = scattered_domain( G1, aj_1_r1, k, X, Y, N1);
% figure(); pcolor(X, Y, real(us1_1)); shading interp; colorbar; title('scattered field no-rereflections, \Gamma_{1} ')
% toc
% 
% UT_2_0 = U_i - us2_0;
% UT_1_1 = U_i - us2_0 - us1_1;
% 
% figure(); pcolor(X, Y, real(UT_2_0)); shading interp; colorbar; title('total field no-rereflections, \Gamma_{2} first ')
% 
% figure(); pcolor(X, Y, real(UT_1_1)); shading interp; colorbar; title('total field 1 re-reflection \Gamma_{2} first')
% 
% %%% Does not work
% %crazy idea:
% UT_12_r1 = U_i - us1_1 - us2_1;
% 
% figure(); pcolor(X, Y, real(UT_12_r1)); shading interp; colorbar; title('Test case, total field 1 re-reflection on each screen')
