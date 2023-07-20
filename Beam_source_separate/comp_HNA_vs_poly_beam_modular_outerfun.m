% In this function we will be looking at the effect of the beam source only
% and making a few adaptations to some of the current general functions.
% This will be largely similar to the modular_outer_fun.m it will change
% when we get to Step 1.


clear all
clear classes

% addpath('/Users/ojp18/Dropbox/Mac/Documents/GitHub/HNA_BEM_Multiple_Scatterers/General_functions')
addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/Multiple_scattering_problems')
addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/General_functions')
addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/BEAM_HNABEMLAB')
addPathsHNA
vertices1 = [-2*pi 2*pi;
    0, 0];
Gamma1=Screen(vertices1);

vertices2 = [2*pi 0;
    5*pi 3*pi];
            
Gamma2=Screen(vertices2);

% General set up
%wavenumber
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

N_approx = 2^(-10);
[x1, y1, t1, t1_mid, h1, h1vector, N1, L1] = discretisation_variables(G1, N_approx, kwave);
[x2, y2, t2, t2_mid, h2, h2vector, N2, L2] = discretisation_variables(G2, N_approx, kwave);

% Infor needed for plotting
x1_plot_1D = [t1_mid(1): 0.001: t1_mid(end)];
x1_plot = linspace(x1_col(1), x1_col(end), length(x1_plot_1D)).';
y1_plot = linspace(y1_col(1), y1_col(end), length(x1_plot_1D)).'; 

x2_plot_1D = [t2_mid(1): 0.001: t2_mid(end)];

% x2_plot_1D = linspace(t2_mid(1), t2_mid(1), 1e4).';
x2_plot = linspace(x2_col(1), x2_col(end), length(x2_plot_1D)).';
y2_plot = linspace(y2_col(1), y2_col(end), length(x2_plot_1D)).';

%%
% Step 0
[v_N_G1_r0] = multiple_Scattering_2screen_step0(kwave, theta, ...
    d, vertices1, L1, col_points1, x1_col, y1_col, h1, t1_mid, t1, colMatrix1, VHNA1, C1, C2);

phi1_0 = @(x) v_N_G1_r0.eval(x, 1) + 2*duidn(vertices1, L1, kwave, d, x);

% phi1_0 = v_N_G1_r0.eval(x1_plot_1D.', 1) + 2*duidn(vertices1, L1, kwave, d, x1_plot_1D.');
figure()
plot(x1_plot_1D/L1, real(phi1_0(x1_plot_1D.')), 'DisplayName', 'OP code', 'LineStyle', '--')
hold on
plot(x1_plot_1D/L1, real(v_N1.eval(x1_plot_1D.',1) + GOA1.eval(x1_plot_1D.', 1)), 'DisplayName', 'AG code', 'LineStyle', '-.')
legend show
xlim([-0.05 1.05])
title('Soluiton on the boundary for r=0, screens in line, long way away')

%% Step 1
N_min = 6;
N_max = 6;
for j = N_min:N_max
    disp(j)
    N_approx = 2^(-j);
    [x1, y1, t1, t1_mid, h1, h1vector, N1, L1] = discretisation_variables(G1, N_approx, kwave);
    [x2, y2, t2, t2_mid, h2, h2vector, N2, L2] = discretisation_variables(G2, N_approx, kwave);


    N_approx_inner = N_approx/2;
    [y1nq_1_inner, y2nq_1_inner, ~, t1_mid_inner, h1_inner, ~, ~, ~] =  discretisation_variables(G1, N_approx_inner, kwave);

% This all needs to be wrapped up in an outer function ideally
% in this case ell = 2, j = 1
    phi_1_outer = phi1_0(t1_mid.');
    phi_1_inner = phi1_0(t1_mid_inner.');
    
    [f_2_1(j, :), f_2_1_uinc(j, :), f_2_1_beam(j, :),...
        ~,S21_phi1_0(j, :), LoB_uinc(j, :), ...
        LoB_beam(j, :),~, ~] ...
    = beam_compute_RHS_vec_given_coll_vec(vertices2, L2, kwave, d, ...
        theta, n1, col_points2, x2_col, y2_col, C1, C2, x1, y1, ...
        h1, phi_1_outer, x2, y2, t2_mid, ...
        t2, h2, y1nq_1_inner, y2nq_1_inner, h1_inner, ...
        phi_1_inner);
end

%% Looking at the beam part of the solution
figure()
for j = N_min:N_max
    plot(col_points2/L2, real(f_2_1_beam(j, :)))
    hold on
    
end
title('Beam source, RHs vec')

figure()
for j = N_min:N_max
    plot(col_points2/L2, real(f_2_1(j, :)))
    hold on
    
end
title('RHS vec combined')

figure()
for j = N_min:N_max
    plot(col_points2/L2, real(S21_phi1_0(j, :)))
    hold on
    
end
title('S21 phi1_0')


%%

v_N_G2_r1 = compute_coeffs_given_A_and_f(colMatrix2, f_2_1(end, :).', VHNA2);

phi_2_r1 = get_phi_j_r(v_N_G2_r1, vertices2, L2, kwave, d, h1, x1, ...
    y1, n2, x2_plot_1D, x2_plot, y2_plot, phi1_0(t1_mid.'));

figure()
plot(x2_plot_1D/L2, real( phi_2_r1 ))
title('Approximation of $\phi_{2}^{(1)}$')

%%

f_test = f_2_1_uinc - f_2_1_beam;
figure()
for j = N_min:N_max
    plot(col_points2/L2, real(f_test(j, :)))
    hold on
    
end
title('RHS vec combined test')

%% computing beam indicent in the field

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
%% 0th step

N_approx = 2^-6;

[x1, y1, t1, t1_mid, h1, h1vector, N1, L1] = discretisation_variables(G1, N_approx, kwave);


phi1_0_vec = phi1_0(t1_mid.');
[us_1_0] = compute_scattered_field_beam(kwave, X, Y, x1, y1, h1, phi1_0_vec);

figure()
pcolor(X, Y, real(us_1_0))
shading interp; colorbar

%% Now looking at beam only problem, beam source is us_1_0, solving only for beam incident not PW.
% v_N_G2_r1_beam_only = compute_coeffs_given_A_and_f(colMatrix2, f_2_1_beam(end, :).', VHNA2);

% LoB = 


%% Complete problem
[x2, y2, t2, t2_mid, h2, h2vector, N2, L2] = discretisation_variables(G2, N_approx, kwave);

phi_2_r1 = get_phi_j_r(v_N_G2_r1, vertices2, L2, kwave, d, h1, x1, ...
    y1, n2, t2_mid, x2, y2, phi1_0(t1_mid.'));

[us_2_1] = compute_scattered_field_beam(kwave, X, Y, x2, y2, h2, phi_2_r1);

figure()
pcolor(X, Y, real(us_2_1))
shading interp; colorbar


figure(); pcolor(X, Y, real(ui - us_1_0)); shading interp; colorbar
figure(); pcolor(X, Y, real(ui - us_1_0 - us_2_1)); shading interp; colorbar

%% Step 2
N_approx = 2^(-6);
[x1, y1, t1, t1_mid, h1, h1vector, N1, L1] = discretisation_variables(G1, N_approx, kwave);
[x2, y2, t2, t2_mid, h2, h2vector, N2, L2] = discretisation_variables(G2, N_approx, kwave);


N_approx_inner = N_approx/2;
[y1nq_2_inner, y2nq_2_inner, ~, t2_mid_inner, h2_inner, ~, ~, ~] =  discretisation_variables(G2, N_approx_inner, kwave);


phi_2_1_outer = get_phi_j_r(v_N_G2_r1, vertices2, L2, kwave, d, h1, x1, ...
y1, n2, t2_mid, x2, y2, phi1_0(t1_mid.'));

phi_2_1_inner = get_phi_j_r(v_N_G2_r1, vertices2, L2, kwave, d, h1, x1, ...
y1, n2, t2_mid_inner, y1nq_2_inner, y2nq_2_inner, phi1_0(t1_mid.'));

[f_1_2, f_1_2_uinc, f_1_2_beam,...
    ~,S12_phi2_1, LoB_uinc_2, ...
    LoB_beam_2,~, ~] ...
= beam_compute_RHS_vec_given_coll_vec(vertices1, L1, kwave, d, ...
    theta, n2, col_points1, x1_col, y1_col, C1, C2, x2, y2, ...
    h2, phi_2_1_outer, x1, y1, t1_mid, ...
    t1, h1, y1nq_2_inner, y2nq_2_inner, h2_inner, ...
    phi_2_1_inner);

%% 
v_N_G1_r2 = compute_coeffs_given_A_and_f(colMatrix1, f_1_2, VHNA1);

phi_1_r2 = get_phi_j_r(v_N_G1_r2, vertices1, L1, kwave, d, h2, x2, ...
    y2, n1, x1_plot_1D, x1_plot, y1_plot, phi_2_1_outer);
%%

N_approx = 2^(-6);
[x1, y1, t1, t1_mid, h1, h1vector, N1, L1] = discretisation_variables(G1, N_approx, kwave);
[x2, y2, t2, t2_mid, h2, h2vector, N2, L2] = discretisation_variables(G2, N_approx, kwave);

figure()
plot(x1_plot_1D/L1, real( phi_1_r2 ))
title('Approximation of $\phi_{1}^{(2)}$')

phi_1_r2 = get_phi_j_r(v_N_G1_r2, vertices1, L1, kwave, d, h2, x2, ...
    y2, n1, t1_mid, x1, y1, phi_2_1_outer);
    
[us_1_2] = compute_scattered_field_beam(kwave, X, Y, x1, y1, h1, phi_1_r2);

figure()
pcolor(X, Y, real(us_2_1))
shading interp; colorbar  


figure(); pcolor(X, Y, real(ui - us_1_0 - us_2_1 - us_1_2)); shading interp; colorbar

%% Step 3 computing phi2_3
phi_1_2_outer = get_phi_j_r(v_N_G1_r2, vertices1, L1, kwave, d, h2, x2, ...
y2, n1, t1_mid, x1, y1, phi_2_1_outer);

phi_2_1_inner = get_phi_j_r(v_N_G2_r1, vertices2, L2, kwave, d, h1, x1, ...
y1, n2, t2_mid, x2, y2, phi1_0(t1_mid.'));

phi_1_2_inner = get_phi_j_r(v_N_G1_r2, vertices1, L1, kwave, d, h2, x2, ...
y2, n1, t1_mid_inner, y1nq_1_inner, y2nq_1_inner, phi_2_1_inner);

[f_2_3, f_2_3_uinc, f_2_3_beam,...
    ~,S21_phi1_2, LoB_uinc_3, ...
    LoB_beam_3,~, ~] ...
= beam_compute_RHS_vec_given_coll_vec(vertices2, L2, kwave, d, ...
    theta, n1, col_points2, x2_col, y2_col, C1, C2, x1, y1, ...
    h1, phi_1_2_outer, x2, y2, t2_mid, ...
    t2, h2, y1nq_1_inner, y2nq_1_inner, h1_inner, ...
    phi_1_2_inner);

v_N_G2_r3 = compute_coeffs_given_A_and_f(colMatrix2, f_2_3, VHNA2);

phi_2_r3 = get_phi_j_r(v_N_G2_r3, vertices2, L2, kwave, d, h1, x1, ...
    y1, n2, x2_plot_1D, x2_plot, y2_plot, phi_1_2_outer);

figure()
plot(x2_plot_1D/L2, real( phi_2_r3 ))
title('Approximation of $\phi_{2}^{(3)}$')

phi_2_r3 = get_phi_j_r(v_N_G2_r3, vertices2, L2, kwave, d, h1, x1, ...
    y1, n2, t2_mid, x2, y2, phi_1_2_outer);
    
[us_2_3] = compute_scattered_field_beam(kwave, X, Y, x2, y2, h2, phi_2_r3);

figure()
pcolor(X, Y, real(us_2_3))
shading interp; colorbar  


% figure(); pcolor(X, Y, real(ui - us_1_0 - us_2_1 - us_1_2 - us_2_3)); shading interp; colorbar



%% various attempts of plotting
% Following is good behind Gamma_{1} not behind Gamma_2.
figure(); pcolor(X, Y, real(ui - us_1_0./2 - us_2_1./2 - us_1_2./2)); shading interp; colorbar


% Following is good for Gamma 2 but not gamma1. Think we need the wave to
% leave the system for it to look reasonable or somehow compute both??

figure(); pcolor(X, Y, real(ui - us_1_0./1 - us_2_1./1 )); shading interp; colorbar

figure(); pcolor(X, Y, real(ui - us_2_1 - us_1_2)); shading interp; colorbar

figure(); pcolor(X, Y, real(ui - us_1_2 - us_2_3)); shading interp; colorbar