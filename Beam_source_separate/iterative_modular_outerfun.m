% outer function computing R iterations

clear all
clear classes

% adding paths, dependent on computer and matlab version
addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/Multiple_scattering_problems')
addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/General_functions')
addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/BEAM_HNABEMLAB')
% addpath('../General_functions/')
% addpath('../Multiple_scattering_problems/')
% addpath('../../BEAM_HNABEMLAB/')
addPathsHNA

R = 20; 
% general geometrical set up
% Case 1
% vertices1 = [-2*pi 2*pi;
%     0, 0];
% 
% vertices2 = [2*pi 0;
%     5*pi 3*pi];

% % Case 2
% vertices1 = [-2*pi 0;
%     0, 0];
% 
% vertices2 = [2*pi 0;
%     3*pi 3*pi];

% Case 3
vertices1 = [-2*pi 2*pi;
    0, 0];

vertices2 = [2*pi 0;
    2*pi 3*pi];

% Case 3a
vertices1 = [-2*pi 2*pi;
    0, 0];

vertices2 = [2*pi 0;
    (2*pi + 0.5) 3*pi];

% Case 4 - Screens facing away from each other - makes a difference which
% order things are put in. How can we fix this?
% vertices1 = [-2*pi -2*pi;
%     0 0];
% 
% vertices2 = [2*pi 0;
%     4*pi -3*pi];

Gamma1=Screen(vertices1);
Gamma2=Screen(vertices2);


kwave=10;
theta = 0;
% theta = pi/4; %% need to change so it is consistent
d = [sin(theta) -cos(theta) ];
uinc=planeWave(kwave,d);

pMax = 4; %polynomial degree
cL = 2; %layers of grading per polynomial degree
sigmaGrad=0.15; %grading ratio
nLayers = cL*(pMax+1)-1; %number of layers of grading
throwAwayParam = 0; %no need to remove any basis elements
OverSample = 1.25; %choose amount to oversample by (40% here)

[v_N1, GOA1, colMatrix1, colRHS1, col_points1,...
v_N2, GOA2, colMatrix2, colRHS2, col_points2, VHNA1, VHNA2] ...
    = AG_code_pulling_out_info(pMax, cL, sigmaGrad, nLayers, OverSample, ...
    Gamma1, Gamma2, kwave, uinc );

% Info needed for our solve
C1 = 1;
C2 = pi;
% Isolating the collocation points in cartesian coordinates
L1 = sqrt( (vertices1(2, 1) - vertices1(1, 1))^2 + (vertices1(2, 2) - vertices1(1, 2))^2 );
x1_col = vertices1(1, 1) + col_points1*(vertices1(2, 1) - vertices1(1, 1))/L1;
y1_col = vertices1(1, 2) + col_points1*( vertices1(2, 2) - vertices1(1, 2) )/L1;

L2 = sqrt( (vertices2(2, 1) - vertices2(1, 1))^2 + (vertices2(2, 2) - vertices2(1, 2))^2 );
x2_col = vertices2(1, 1) + col_points2*(vertices2(2, 1) - vertices2(1, 1))/L2;
y2_col = vertices2(1, 2) + col_points2*( vertices2(2, 2) - vertices2(1, 2) )/L2;

G1 = [vertices1(1, 1), vertices1(1, 2), vertices1(2, 1), vertices1(2, 2)];
G2 = [vertices2(1, 1), vertices2(1, 2), vertices2(2, 1), vertices2(2, 2)];

n1 = [-(G1(4) - G1(2)), G1(3) - G1(1)]/L1;
n2 = [-(G2(4) - G2(2)), G2(3) - G2(1)]/L2;

[dot(n1, d), dot(n2, d)]

keyboard

N_approx = 2^(-6);

[aj_1_r, aj_2_r, phi1_r_outer, phi2_r_outer] = ...
    compute_coeff_LOB_for_R_iterations(kwave, N_approx, G1, G2, ...
    vertices1, vertices2, R, theta,...
    col_points1, x1_col, y1_col, col_points2, x2_col, y2_col, VHNA1,...
    VHNA2, colMatrix1, colMatrix2, d, n1, n2, C1, C2);

phi1_x = linspace(0, 1, length(phi1_r_outer(1, :)));
phi2_x = linspace(0, 1, length(phi2_r_outer(1, :)));

%%
R_min_plot = 1;
R_max_plot = 4;
figure()
for r = R_min_plot:R_max_plot
    subplot(2, 1, 1)
    plot(phi1_x, real(phi1_r_outer(r, :)), 'DisplayName', strcat('r = ', num2str(2*r - 2)))
    hold on
    
    subplot(2, 1, 2)
    plot(phi2_x, real(phi2_r_outer(r, :)), 'DisplayName', strcat('r = ', num2str(2*r - 1)))
    hold on
    
end

subplot(2, 1, 1)
title('Approximation of $\phi_{1}^{(2r)}$', 'fontsize',18,'interpreter','latex')
xlabel('$s/L$', 'fontsize',18,'interpreter','latex')
ylabel('$\phi_{1}^{(2r)}(s)$', 'fontsize',18,'interpreter','latex')
xlim([-0.05 1.05])
legend show

subplot(2, 1, 2)
title('Approximation of $\phi_{2}^{(2r+1)}$', 'fontsize',18,'interpreter','latex')
xlabel('$s/L$', 'fontsize',18,'interpreter','latex')
ylabel('$\phi_{2}^{(2r+1)}(s)$', 'fontsize',18,'interpreter','latex')
xlim([-0.05, 1.05])
legend show

%% Plotting in the domain
[x1, y1, t1, t1_mid, h1, h1vector, N1, L1] = discretisation_variables(G1, N_approx, kwave);
[x2, y2, t2, t2_mid, h2, h2vector, N2, L2] = discretisation_variables(G2, N_approx, kwave);

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

%% Computing incident plane wave
ui = incident2d(kwave, theta, X, Y);

% figure()
% pcolor(X, Y, real(ui))
% shading interp; colorbar

[us_phi1, x1, y1] = ...
    compute_scattered_field_beam(kwave, X, Y, x1, y1, h1,...
    phi1_r_outer(end, :).', G1, t1_mid);

[us_phi2, x2, y2] = ...
    compute_scattered_field_beam(kwave, X, Y, x2, y2, h2, ...
    phi2_r_outer(end, :).', G2, t2_mid);

u_best_approx = ui - (us_phi1 + us_phi2);

figure()
pcolor(X, Y, real(u_best_approx));
hold on
Gamma_1 = plot(x1, y1);
Gamma_1.LineWidth = 2;
Gamma_1.Color = [0 0 0];
Gamma_2 = plot(x2, y2);
Gamma_2.LineWidth = 2;
Gamma_2.Color = [0 0 0];
shading interp; colorbar
txt_title = ['Total solution in the field, $k = $', num2str(kwave),...
    '\theta = ', num2str(theta), '$R = $', num2str(R)];
title('Total solution in the field, $k = 10$, $\theta = 0$, R = 20.',...
    'fontsize',18)
xlabel('$x$', 'fontsize',18)
ylabel('$y$', 'fontsize',18)