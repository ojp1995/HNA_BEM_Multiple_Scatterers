% Here we will be comparing the polynomial code with the HNA code

clear all

load('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/Poly_approx_space_PIM/polycode_S21_phi1_0.mat')

addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/general_tests')
addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/Multiple_scattering_problems')
addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/General_functions')
addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/BEAM_HNABEMLAB/')
addPathsHNA

S21phi1_0_poly = S21_phi1_0;
phi1_0_poly = phi1_0;

vertices1 = [G1(1), G1(2) ;
    G1(3), G1(4)];
Gamma1 = Screen(vertices1);

vertices2 = [G2(1), G2(2) ;
    G2(3), G2(4)];
Gamma2 = Screen(vertices2);

kwave = k;

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

N_approx = C_wl;

[x1, y1, t1, t1_mid, h1, h1vector, N1, L1] = discretisation_variables(G1, N_approx, kwave);
[x2, y2, t2, t2_mid, h2, h2vector, N2, L2] = discretisation_variables(G2, N_approx, kwave);


% Infor needed for plotting
x1_plot_1D = [t1_mid(1): 0.001: t1_mid(end)];
x1_plot = linspace(x1_col(1), x1_col(end), length(x1_plot_1D)).';
y1_plot = linspace(y1_col(1), y1_col(end), length(x1_plot_1D)).'; 

x2_plot_1D = [t2_mid(1): 0.001: t2_mid(end)];

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


% x2_plot_1D = linspace(t2_mid(1), t2_mid(1), 1e4).';
x2_plot = linspace(x2_col(1), x2_col(end), length(x2_plot_1D)).';
y2_plot = linspace(y2_col(1), y2_col(end), length(x2_plot_1D)).';

N_approx_inner = N_approx/2;
[y1nq_1_inner, y2nq_1_inner, ~, t1_mid_inner, h1_inner, ~, ~, ~] =  discretisation_variables(G1, N_approx_inner, kwave);

% This all needs to be wrapped up in an outer function ideally
% in this case ell = 2, j = 1
phi_1_outer = phi1_0(t1_mid.');
phi_1_inner = phi1_0(t1_mid_inner.');
[f_2_r1, ~, S21phi1_0, S22Psi_2_1]...
= compute_RHS_vec_given_coll_vec(vertices2, L2, kwave, d, ...
    theta, n1, col_points2, x2_col, y2_col, C1, C2, x1, y1, ...
    h1, phi_1_outer, x2, y2, t2_mid, ...
    t2, h2, y1nq_1_inner, y2nq_1_inner, h1_inner, ...
    phi_1_inner);

v_N_G2_r1 = compute_coeffs_given_A_and_f(colMatrix2, f_2_r1, VHNA2);

phi_2_r1 = get_phi_j_r(v_N_G2_r1, vertices2, L2, kwave, d, h1, x1, ...
    y1, n2, x2_plot_1D, x2_plot, y2_plot, phi1_0(t1_mid.'));

figure()
plot(x2_plot_1D/L2, real( phi_2_r1 ))
title('Approximation of $\phi_{2}^{(1)}$')

%% comparison between two methods

[~, ~, S21phi1_0_HNA, ~, ~] = compute_RHS_vec_given_coll_vec(vertices2, L2, kwave, d, ...
    theta, n1, t2_mid, x2, y2, C1, C2, x1, y1, ...
    h1, phi_1_outer, x2, y2, t2_mid, ...
    t2, h2, y1nq_1_inner, y2nq_1_inner, h1_inner, ...
    phi_1_inner);


err_poly_HNA_S21phi1_0 = abs(S21phi1_0_poly - S21phi1_0_HNA);

figure()
plot(t2_mid/L2, real(S21phi1_0_HNA), 'DisplayName', 'HNA code')
hold on
plot(t2_mid/L2, real(S21phi1_0_poly), 'DisplayName', 'Poly code')
hold off
legend show
title('Approximation of $S_{21}\phi_{1}^{(0)}$ using different approximation methods')
xlabel('s/L')
ylabel('Approximation to $S_{21}\phi_{1}^{(0)}(s)$')


figure()
plot(t2_mid/L2, err_poly_HNA_S21phi1_0)
title('Difference betweeen the polynomial approximation and HNA')
xlabel('s/L')
ylabel('Difference')

%% Comparing the approximatiosn to phi1_0
err_phi1_0_poly_vs_HNA = abs(phi1_0_poly - phi_1_outer);

figure()
plot(t1_mid/L1, err_phi1_0_poly_vs_HNA)

figure()
plot(t1_mid(5:end - 5)/L1, err_phi1_0_poly_vs_HNA(5:end - 5))

%% computing S21phi1_0 using poly S21 and HNA phi1_0

S21poly_phi1_0_HNA = S21*phi_1_outer;

err_poly_vs_HNA_S21phi1_0_adapted = abs(S21poly_phi1_0_HNA - S21phi1_0_HNA);

figure()
plot(t2_mid/L2, err_poly_vs_HNA_S21phi1_0_adapted)
xlabel('$s/L_{2}$', 'fontsize',15,'interpreter','latex')
ylabel('Absolute difference', 'fontsize',15,'interpreter','latex')
title('Absolute difference between the HNA code and polynomial code of $ \frac{i}{4} \int_{0}^{L_{1}} H_{0}^{(1)}(k \vert \mathbf{x}_{1}(s) - \mathbf{y}(t) ) \phi_{1}^{(0)}(t) \mathrm{d} t$', 'fontsize',15,'interpreter','latex')
xlim([-0.1 1.1])

%% Max compairison

[max(real(S21phi1_0_HNA)) max(real(S21poly_phi1_0_HNA)) max(real(S21phi1_0_poly))]
