% Outer function to test new more modular approach of solving the problem

clear all
clear classes

% addpath('/Users/ojp18/Dropbox/Mac/Documents/GitHub/HNA_BEM_Multiple_Scatterers/General_functions')
addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/BEAM_HNABEMLAB')
addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/General_functions')
addPathsHNA
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

%%
% % Step 1 - Convergence test
N_min = 1;
N_max = 8;
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
    [f_2_r1(j, :), ~, S21phi1_0(j, :), S22Psi_2_1(j, :)]...
    = compute_RHS_vec_given_coll_vec(vertices2, L2, kwave, d, ...
        theta, n1, col_points2, x2_col, y2_col, C1, C2, x1, y1, ...
        h1, phi_1_outer, x2, y2, t2_mid, ...
        t2, h2, y1nq_1_inner, y2nq_1_inner, h1_inner, ...
        phi_1_inner);

end

%%
% this part of the function is given f and A compute the coefficients and
% then ideally provide a function for approximating 
v_N_G2_r1 = compute_coeffs_given_A_and_f(colMatrix2, f_2_r1(end, :).', VHNA2);

% phi2_r1 = @(t2_1da, x2a, y2a, t1_1da) v_N_G2_r1.eval(t2_1da, 1) ...
%     + 2*duidn(vertices2, L2, kwave, d, t2_1da).'...
%     + midpoint_dphikdn_f_diff_screen(kwave, x2a, y2a, h1, x1, ...
%     y1, phi1_0(t1_1da.').', n2);

phi_2_r1 = get_phi_j_r(v_N_G2_r1, vertices2, L2, kwave, d, h1, x1, ...
    y1, n2, x2_plot_1D, x2_plot, y2_plot, phi1_0(t1_mid.'));

figure()
plot(x2_plot_1D/L2, real( phi_2_r1 ))
title('Approximation of $\phi_{2}^{(1)}$')

%%
% error computation
% we are wanting the overall error and the max values in each
figure()
N_min_plot = 5;

for j = N_min:N_max
    err_f_2_r1(j) = sum((abs(f_2_r1(end, :) - f_2_r1(j, :)))./f_2_r1(end, :))/length(col_points2);
%     err_ddn_S21phi1_0(j) = sum((ddn_S21phi1_0(end, :) - ddn_S21phi1_0(j, :))./ddn_S21phi1_0(end, :));
    err_S21phi1_0(j) = sum((abs(S21phi1_0(end, :) - S21phi1_0(j, :))./S21phi1_0(end, :)))/length(col_points2);
    err_S22Psi_2_1(j) = sum((abs(S22Psi_2_1(end, :) - S22Psi_2_1(j, :))./S22Psi_2_1(end, :)))/length(col_points2);

    % max values real
    max_real_f_2_r1(j) = max(real(f_2_r1(j, :)));
%     max__real_ddn_S21phi1_0(j) = max(real(ddn_S21phi1_0(j, :)));
    max_real_S21phi1_0(j) = max(real(S21phi1_0(j, :)));
    max_real_S22Psi_2_1(j) = max(real(S22Psi_2_1(j, :)));
    
   
 
end

for j = N_min_plot:N_max
     txt1 = [ 'dof = 2^', num2str(j)];
    plot(col_points2/L2, real(f_2_r1(j, :)), 'DisplayName', txt1)
    hold on

end
title('RHS vector $f_{2}^{(1)}$')
hold off
legend show

figure()
for l = N_min_plot:N_max
    txt2 = ['dof = 2^', num2str(l)];
    plot(col_points2/L2, real(S21phi1_0(l, :)), 'DisplayName', txt2)
    hold on
    
end
title('Plot of $S_{21} \phi_{1}^{(0)}$ for increasing order dof')
hold off
legend show


figure()
for xx= N_min_plot:N_max

    txt3 = ['dof = 2^', num2str(xx)];
    plot(col_points2/L2, real(S22Psi_2_1(xx, :)), 'DisplayName', txt3)
    hold on

end
title('Plot of $S_{22} \Psi_{2}^{(1)}$ for increasing dof')
hold off
legend show

% estimated order of convergence

for j = N_min+1:N_max - 1

    EOC_f_2_r1(j) = log2(err_f_2_r1(j-1)/err_f_2_r1(j));

    EOC_S21phi1_0(j) = log2(err_S21phi1_0(j-1)/err_S21phi1_0(j));

    EOC_S22Psi_2_1(j) = log2(err_S22Psi_2_1(j-1)/err_S22Psi_2_1(j));

end

[err_f_2_r1.', err_S21phi1_0.', err_S22Psi_2_1.']

[EOC_f_2_r1.', EOC_S21phi1_0.', EOC_S22Psi_2_1.']

[max_real_f_2_r1.', max_real_S21phi1_0.', max_real_S22Psi_2_1.']




















