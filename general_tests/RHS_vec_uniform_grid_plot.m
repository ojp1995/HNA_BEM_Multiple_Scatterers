% testing to see the behaviour of the right hand side vector, f_2^(1) at a
% larger number of points to see more of its behaviour.
%
% We will also try to compare this to complete polynomial approximation
% space to see how closely it matches.

clear all

% first we need to creata  generic set up to then plot the righ hand sides:

clear all
clear classes

addpath('/Users/ojp18/Dropbox/Mac/Documents/GitHub/HNA_BEM_Multiple_Scatterers/General_functions')
addpath('../Multiple_scattering_problems/')

% TODO: check what order the coefficients need to be in!
% General set up

% Geometric set up
% vertices1 = [-2*pi 0;
%     0 0];
% 
vertices1 = [-2*pi 2*pi;
    0, 0];
Gamma1=Screen(vertices1);


% vertices2 = [200*pi 0;
%     202*pi 0];

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

disp(length(col_points1))

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

% dui_dn1 = @(nodes) duidn(vertices1, L1, kwave, d, nodes);
% dui_dn2 = @(nodes) duidn(vertices2, L2, kwave, d, nodes);

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
N_f1_test = 1000;
h_f1_test = L1/N_f1_test;
col_points1_test = [0:h_f1_test:L1];
x1_col_test = vertices1(1, 1) + col_points1_test*(vertices1(2, 1) - vertices1(1, 1))/L1;
y1_col_test = vertices1(1, 2) + col_points1_test*( vertices1(2, 2) - vertices1(1, 2) )/L1;


N_step0_test = 11;
figure()
for j = N_step0_test:N_step0_test
    N_approx = 2^(-j);
    [x1, y1, t1, t1_mid, h1, h1vector, N1, L1] = discretisation_variables(G1, N_approx, kwave);
%     [x2, y2, t2, t2_mid, h2, h2vector, N2, L2] = discretisation_variables(G2, N_approx, kwave);
    f1_0(:, 1) = Step0_compute_RHS_given_coll_vec(kwave, theta, ...
    col_points1_test, x1_col_test, y1_col_test, L1, vertices1, d, h1, t1_mid, t1, C1, C2 );

    
%     [v_N_G1_r0, f1_0] = multiple_Scattering_2screen_step0(kwave, theta, ...
%         d, vertices1, L1, col_points1, x1_col, y1_col, h1, t1_mid, t1, colMatrix1, VHNA1, C1, C2);
    txt = ['#dof = 2^-', num2str(j)];
    plot(col_points1_test, real(f1_0), 'DisplayName', txt)
    hold on
    
end
title('Increasing quadrature approx plot of RHS vector $f_{1}^{(0)}$. Fixed collocation points')
legend show

%%
% v_N_G1_r0 = Step0_compute_coeffs_given_A_and_f(colMatrix1, f1_0, VHNA1);
phi1_0 = @(x) v_N_G1_r0.eval(x, 1) + 2*duidn(vertices1, L1, kwave, d, x);

% phi1_0 = v_N_G1_r0.eval(x1_plot_1D.', 1) + 2*duidn(vertices1, L1, kwave, d, x1_plot_1D.');
figure()
plot(x1_plot_1D/L1, real(phi1_0(x1_plot_1D.')), 'DisplayName', 'OP code', 'LineStyle', '--')
hold on
plot(x1_plot_1D/L1, real(v_N1.eval(x1_plot_1D.',1) + GOA1.eval(x1_plot_1D.', 1)), 'DisplayName', 'AG code', 'LineStyle', '-.')
legend show
xlim([-0.05 1.05])
title('Soluiton on the boundary for r=0, screens in line, long way away')

[x2, y2, t2, t2_mid, h2, h2vector, N2, L2] = discretisation_variables(G2, N_approx, kwave);


%% computing the right hand side vector at more points that normal
% So we are going to edit the inputs s_ell, x1s_ell, y1s_ell, in the 
% function compute_RHS_vec_given_coll_vec.m so that we can plot f_2^(1) at
% more points and other parts of it.
% 
% To start with we will use a uniform grid.
keyboard
N_RHS = 1000;

h_RHS = L2/N_RHS;
s_RHS = [0: h_RHS: L2 ];
x2_col_RHS = vertices2(1, 1) + s_RHS*(vertices2(2, 1) - vertices2(1, 1))/L2;
y2_col_RHS = vertices2(1, 2) + s_RHS*( vertices2(2, 2) - vertices2(1, 2) )/L2;

N_approx_inner = N_approx/2;
[y1nq_1_inner, y2nq_1_inner, ~, t1_mid_inner, h1_inner, ~, ~, ~] =  discretisation_variables(G1, N_approx_inner, kwave);

phi_1_outer = phi1_0(t1_mid.');
phi_1_inner = phi1_0(t1_mid_inner.');

[f_2_r1, ddn_S21phi1_0, S21phi1_0, S22Psi_2_1, Psi21]...
    = compute_RHS_vec_given_coll_vec(vertices2, L2, kwave, d, ...
        theta, n1, s_RHS, x2_col_RHS, y2_col_RHS, C1, C2, x1, y1, ...
        h1, phi_1_outer, x2, y2, t2_mid, ...
        t2, h2, y1nq_1_inner, y2nq_1_inner, h1_inner, ...
        phi_1_inner);

figure()
plot(s_RHS/L2, real(f_2_r1))
title('Plot of $f_{2}^{(1)}$ on a uniform grid of points')

figure()
plot(s_RHS/L2, real(S21phi1_0))
title('Plot of $S_{21}\phi_{1}^{(0)}$ on a uniform grid of points')

figure()
plot(s_RHS/L2, real(S22Psi_2_1))
title('Plot of $S_{22} \Psi_{2}^{(1)}$ on a uniform grid of points')

figure()
plot(t2_mid/L2, real(ddn_S21phi1_0))
title('Plot of $\frac{\partial}{\partial n} S_{22}\phi_{1}^{0}$')

figure()
plot(t2_mid/L2, real(Psi21))
title('Plot of $\Psi_{2}^{(1)}$')

