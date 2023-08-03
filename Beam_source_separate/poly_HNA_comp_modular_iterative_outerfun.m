% outerfunction computing R iterations.

clear all
clear classes

% load('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/Poly_approx_space_PIM/polycode_test1_R_20_dof2_6.mat')
load('../Poly_approx_space_PIM/polycode_test_knifeedge_R_20_dof2_6.mat')
phi1_r_poly = phi1_r;
phi2_r_poly = phi2_r;
% addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/Multiple_scattering_problems')
% addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/General_functions')
% addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/general_tests')
% addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/BEAM_HNABEMLAB')
addpath('../general_tests/')
addpath('../General_functions/')
addpath('../Multiple_scattering_problems/')
addpath('../../BEAM_HNABEMLAB/')
addpath('../Poly_approx_space_PIM/')
addPathsHNA

vertices1 = [G1(1), G1(2) ;
    G1(3), G1(4)];
Gamma1 = Screen(vertices1);

vertices2 = [G2(1), G2(2) ;
    G2(3), G2(4)];
Gamma2 = Screen(vertices2);

kwave = k;

% vertices1 = [-2*pi 2*pi;
%     0, 0];
% Gamma1=Screen(vertices1);

% vertices2 = [2*pi 0;
%     5*pi 3*pi];
% Gamma2=Screen(vertices2);

% General set up
%wavenumber
% kwave=10;
% theta = 0;
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

% R = 20; 
%%
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

%% compare to the polynomial code 

% intergation nodes for L1 error computations
L1_temp = L1;

N_L1_err = 2^12;
h1_L1_err = L1/N_L1_err;
t1_L1_err_grid = [0: h1_L1_err: L1];
tq1 = [h1_L1_err/2 : h1_L1_err: L1 - h1_L1_err];
% computing the midpoints coordinates for \Gamma
xq1 = G1(1) + tq1*(G1(3) - G1(1))/(L1);
yq1 = G1(2) + tq1*(G1(4) - G1(2))/(L1);

% xq1 = G1(1)+((1:N_L1_err)-0.5)*(G1(3)-G1(1))/N_L1_err;
% yq1 = G1(2)+((1:N_L1_err)-0.5)*(G1(4)-G1(2))/N_L1_err;

% inner variables
[x2_in, y2_in, ~, ~, h2_in, ~, ~, ~] = ...
    discretisation_variables(G2, N_approx, kwave);

% N_L1_err = 2^6;
h2_L1_err = L2/N_L1_err;
t2_L1_err_grid = [0: h2_L1_err: L2];
tq2 = [h2_L1_err/2 : h2_L1_err: L2 - h2_L1_err];
% computing the midpoints coordinates for \Gamma
xq2 = G2(1) + tq2*(G2(3) - G2(1))/(L2);
yq2 = G2(2) + tq2*(G2(4) - G2(2))/(L2);

% xq1 = G1(1)+((1:N_L1_err)-0.5)*(G1(3)-G1(1))/N_L1_err;
% yq1 = G1(2)+((1:N_L1_err)-0.5)*(G1(4)-G1(2))/N_L1_err;

% inner variables
[x1_in, y1_in, ~, ~, h1_in, ~, ~, ~] = ...
    discretisation_variables(G1, N_approx, kwave);


% error computations for last iteration
[L1_err_G1_no_norm, L1_err_norm_poly_true_G1, L1_err_norm_HNA_true_G1, phi_poly1, phi_HNA1] = L1_err_poly_HNA(tq1, xq1, yq1, h1_L1_err, vertices1, L1, k, d, n1,...
    phi1_r_poly(end, :), aj_1_r{end},  h2_in, x2_in, ...
    y2_in, phi2_r_outer(end - 1, :).');

[L1_err_G2_no_norm, L1_err_norm_poly_true_G2, L1_err_norm_HNA_true_G2, phi_poly2, phi_HNA2] = L1_err_poly_HNA(tq2, xq2, yq2, h2_L1_err, vertices2, L2, k, d, n2,...
    phi2_r_poly(end, :), aj_2_r{end},  h1_in, x1_in, ...
    y1_in, phi1_r_outer(end - 1, :).');

figure()
plot(tq1/L1, real(phi_poly1), 'DisplayName', 'Poly approx')
hold on
plot(tq1/L1, real(phi_HNA1), 'DisplayName', 'HNA approx')
legend show
xlim([-0.05 1.05])

figure()
plot(tq2/L2, real(phi_poly2), 'DisplayName', 'Poly approx')
hold on
plot(tq2/L2, real(phi_HNA2), 'DisplayName', 'HNA approx')
legend show
xlim([-0.05 1.05])

L1_err_G1_no_norm, L1_err_norm_poly_true_G1, L1_err_norm_HNA_true_G1
L1_err_G2_no_norm, L1_err_norm_poly_true_G2, L1_err_norm_HNA_true_G2
%%

L1 = L1_temp;
% error calculations
for r = 1:R
    
    err_HNAvspoly_phi1(r, :) = abs(phi1_r_outer(r, :) - phi1_r_poly(r, :));
    sum_err_phi1(r) = sum(err_HNAvspoly_phi1(r, :))/length(phi1_r_outer(r, :));
   
    
    err_HNAvspoly_phi2(r, :) = abs(phi2_r_outer(r, :) - phi2_r_poly(r, :));
    sum_err_phi2(r) = sum(err_HNAvspoly_phi2(r, :))/length(phi2_r_outer(r, :));
    
end

% first compare directly with plots one over the other
% comparing \phi1
R_min_plot = 1;
R_max_plot = 3;
figure()
for r = R_min_plot:R_max_plot
    plot(phi1_x, real(phi1_r_outer(r, :)), '-.',  'DisplayName', strcat('HNA, r = ', num2str(2*r - 2)));
    hold on
    plot(phi1_x, real(phi1_r_poly(r, :)), 'DisplayName', strcat('poly, r = ', num2str(2*r - 2)),  'LineStyle', '--');
    

end
plot(phi1_x, real(phi1_r_poly(end, :)), 'DisplayName', strcat('poly, r = ', num2str(2*R - 2)),  'LineStyle', '--');
legend show
title('Comparison between poly and HNA aproximation of $\phi_{1}^{(r)}$', 'fontsize',18,'interpreter','latex')
xlabel('$s/L$', 'fontsize',18,'interpreter','latex')
ylabel('$\phi_{1}^{r}(s)$', 'fontsize',18,'interpreter','latex')
xlim([-0.05 1.05])

figure()
for r = R_min_plot:R_max_plot
    plot(phi2_x, real(phi2_r_outer(r, :)), '-.', 'DisplayName', strcat('HNA, r = ', num2str(2*r - 1)));
    hold on
    plot(phi2_x, real(phi2_r_poly(r, :)), 'DisplayName', strcat('poly, r = ', num2str(2*r - 1)), 'LineStyle', '--');
    
end
plot(phi2_x, real(phi2_r_poly(end, :)), 'DisplayName', strcat('poly, r = ', num2str(2*R - 1)), 'LineStyle', '--');
legend show
title('Comparison between poly and HNA aproximation of $\phi_{2}^{(r)}$', 'fontsize',18,'interpreter','latex')
xlabel('$s/L$', 'fontsize',18,'interpreter','latex')
ylabel('$\phi_{2}^{r}(s)$', 'fontsize',18,'interpreter','latex')
xlim([-0.05 1.05])

% Now compute and plot errors
figure()
for r = R_min_plot:R_max_plot
    plot(phi1_x, real(err_HNAvspoly_phi1(r, :)), 'DisplayName', strcat('r = ', num2str(2*r - 2)));
    hold on
    
end
legend show
title('Absolute error between HNA and polynomial approximation of $\phi_{1}^{(r)}$', 'fontsize',18,'interpreter','latex')
xlabel('$s/L$', 'fontsize',18,'interpreter','latex')
ylabel('Difference', 'fontsize',18,'interpreter','latex')
xlim([-0.05 1.05])

figure()
for r = R_min_plot:R_max_plot
    plot(phi2_x, real(err_HNAvspoly_phi2(r, :)), 'DisplayName', strcat('r = ', num2str(2*r - 1)));
    hold on
    
end
legend show
title('Absolute error between HNA and polynomial approximation of $\phi_{2}^{(r)}$', 'fontsize',18,'interpreter','latex')
xlabel('$s/L$', 'fontsize',18,'interpreter','latex')
ylabel('Difference', 'fontsize',18,'interpreter','latex')
xlim([-0.05 1.05])

figure()
plot(sum_err_phi1, 'DisplayName', 'Averaged error in $\phi_{1}^{(r)}$')
hold on
plot(sum_err_phi2, 'DisplayName', 'Averaged error in $\phi_{2}^{(r)}$')
legend show
title('Summed average between HNA and poly computation for different R', 'fontsize',18,'interpreter','latex')
xlabel('r', 'fontsize',18,'interpreter','latex')
ylabel('Summed difference', 'fontsize',18,'interpreter','latex')

%% Plotting in the domain
error('Why??')
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

figure()
pcolor(X, Y, real(ui))
shading interp; colorbar

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
txt_title = ['Total solution in the field, $k = $', num2str(k),...
    '\theta = ', num2str(theta), '$R = $', num2str(R)];
title('Total solution in the field, $k = 10$, $\theta = 0$, R = 20.',...
    'fontsize',18)
xlabel('$x$', 'fontsize',18)
ylabel('$y$', 'fontsize',18)
%%
error('Seriously slow code ahead!')

us_phi1_r = zeros(length(X), length(Y), R+1);
us_phi2_r = zeros(length(X), length(Y), R+1);
u_r = zeros(length(X), length(Y), R+1);
for r = 1:R+1
    tic
    us_phi1_r(:, :, r) =  ...
        compute_scattered_field_beam(kwave, X, Y, x1, y1, h1, phi1_r_outer(r, :).');
    
    us_phi2_r(:, :, r) =  ...
        compute_scattered_field_beam(kwave, X, Y, x2, y2, h2, phi2_r_outer(r, :).');
    
    u_r(:, :, r) = ui - ( us_phi1_r(:, :, r) + us_phi2_r(:, :, r) );
    toc
    
    figure()
    pcolor(X, Y, real(u_r(:, :, r)));
    shading interp; colorbar
    title(['Total solution in the domain with r = ', num2str(r)], 'fontsize',18,'interpreter','latex')

end

