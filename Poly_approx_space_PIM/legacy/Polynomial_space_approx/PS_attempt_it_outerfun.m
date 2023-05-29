% iterative method using a point source attmept. Method is using a PIM to
% evaluate the singular integrals using a midpoint rule.

clear all

tic

% point source
xps = pi;
yps = 4*pi;

% introducing the screens
G1 = [-2*pi, 2*pi, 0, 0]; %have switched as a test
L1 = sqrt( (G1(3) - G1(1))^2 + (G1(4) - G1(2))^2 );  % length of G1

G2 = [ 2*pi, 0, 5*pi, 3*pi];
L2 = sqrt( (G2(3) - G2(1))^2 + (G2(4) - G2(2))^2 );  % length of G2

C_wl= 1/40

k = 5;  % wavenumber

N1 = ceil(k*L1./(C_wl*2*pi)) % number of itervals on G1
N2 = ceil(k*L2./(C_wl*2*pi)) % number of intervals on G2

theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

% I think the set up has now been introduced so we can start to compute the
% matrices for \Gamma_{1} and \Gamma_{2}

Q = 1;  % quadrature discretisation

[A1, x_col1 A1_sing, A1_smooth] = screen_mat_poly_approx_PIM(N1, 0, L1, k, C1, C2, Q);
[A2, x_col2, A2_sing, A2_smooth] = screen_mat_poly_approx_PIM(N2, 0, L2, k, C1, C2, Q);

[S21, col_points_1] = S21_op(G1, G2, k, N1, N2);
[S12, col_points_2] = S21_op(G2, G1, k, N2, N1);

% automating the re-reflections
R = 4; % this is the number of re-reflections we want
% initialising
u_i_1_r = zeros(length(x_col1), R);
aj_1_r = zeros(length(x_col1), R);
u_i_2_r = zeros(length(x_col2), R);
aj_2_r = zeros(length(x_col2), R);
% phi1_r = zeros(N_sample, R);
% phi2_r = zeros(N_sample, R);


% r = 1 case:
u_i_1_r(:, 1) = PS_incident(G1, L1, k, x_col1, xps, yps);
aj_1_r(:, 1) = A1\u_i_1_r(:, 1);
% [phi1_r(:, 1), x_sample1] = coeff_2soln_midpoint(aj_1_r(:, 1), L1, N_sample-1, N1);

u_i_2_r(:, 1) =  PS_incident(G2, L2, k, x_col2, xps, yps);
aj_2_r(:, 1) = A2\u_i_2_r(:, 1);
% [phi2_r(:, 1), x_sample2] = coeff_2soln_midpoint(aj_2_r(:, 1), L2, N_sample-1, N2);

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
    u_i_1_r(:, r) =   PS_incident(G1, L1, k, x_col1, xps, yps).' - S12*coeff_2soln_midpoint(aj_2_r(:, r-1), L2, N2-1, N2);
%     conpute the coefficients
    aj_1_r(:, r) = A1\u_i_1_r(:, r);
%     compute the solution
%     [phi1_r(:, r), ~] = coeff_2soln_midpoint(aj_1_r(:, r), L1, N_sample-1, N1);
    figure(5)
    plot(x_col1/L1, aj_1_r(:, r), 'DisplayName', sprintf('r = %g', 2*r - 2))
    
%     compute the odd solutions on \Gamma_{2}
%         first compute the incident:
    u_i_2_r(:, r) =  PS_incident(G2, L2, k, x_col2, xps, yps).' - S21*coeff_2soln_midpoint(aj_1_r(:, r), L1, N1-1, N1);
%     conpute the coefficients
    aj_2_r(:, r) = A2\u_i_2_r(:, r);
%     compute the solution
%     [phi2_r(:, r), ~] = coeff_2soln_midpoint(aj_2_r(:, r), L2, N_sample-1, N2);
    
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
ylabel('Re $\left(\phi_{1}^{r} \right)$', 'FontSize', 17, 'interpreter', 'latex')
xlim([-0.1 1.1])
legend show
set(legend,'fontsize',17);
subplot(2, 1, 2)
% title('Solution on \Gamma_{2} for increasing number of iterations')
xlabel('$s/L_{2}$', 'FontSize', 17, 'interpreter', 'latex')
ylabel('Re $ \left(\phi_{2}^{r} \right) $', 'FontSize', 17, 'interpreter', 'latex')
xlim([-0.1 1.1])
legend show
set(legend,'fontsize',17);

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

N_d = 100;

min_leftorbottom = min(X_min, Y_min);
max_toporright = max(X_max, Y_max);
h_d = N_diff/N_d;
% need to change this bit so that it is closer zoomed in on the picture
X = [min_leftorbottom - 5: h_d: max_toporright + 5];
Y = [min_leftorbottom - 5: h_d: max_toporright + 5];
% first compute the incident field

U_i = PS_incident_domain(xps, yps, X, Y, k);

figure(); pcolor(X, Y, real(U_i)); shading interp; title('Point source attempt')

toc

% compute the scattered field
Q_d = 110;

tic
[us1_0, x1, y1] = scattered_domain( G1, aj_1_r(:, 1), k, X, Y, N1);
figure(); pcolor(X, Y, real(us1_0)); shading interp; colorbar; title('scattered field no-rereflections ')
toc

tic
[us2_1, x2, y2] = scattered_domain( G2, aj_2_r(:, 1), k, X, Y, N2);
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

figure(); pcolor(X, Y, real(UT_1_0)); shading interp; colormap(jet); colorbar; title('total field no-rereflections ')
hold on
Gamma_1 = plot(x1, y1);
Gamma_1.LineWidth = 4;
Gamma_1.Color = [0 0 0];

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
