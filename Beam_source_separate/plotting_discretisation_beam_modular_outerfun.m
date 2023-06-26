% outerfunction computing R iterations.  And also plotting the
% discretisation points to investigate the possibility of 

clear all
clear classes

% load('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/Poly_approx_space_PIM/polycode_test1_R_20.mat')

% addpath('/Users/ojp18/Dropbox/Mac/Documents/GitHub/HNA_BEM_Multiple_Scatterers/General_functions')
% addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/Multiple_scattering_problems')
% addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/General_functions')
% addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/BEAM_HNABEMLAB')
addpath('../General_functions/')
addpath('../Multiple_scattering_problems/')
addpath('../../BEAM_HNABEMLAB/')
addPathsHNA

vertices1 = [-2*pi 2*pi;
    0, 0];
Gamma1=Screen(vertices1);

vertices2 = [2*pi 0;
    5*pi 3*pi];
            
Gamma2=Screen(vertices2);

R = 10;

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

%% isolating the discretisation points

mesh1_p = VHNA1.mesh{1}{1}.points(end - 2:end);
mesh1_n = VHNA1.mesh{1}{2}.points(1:3);

mesh2_p = VHNA2.mesh{1}{1}.points(end - 2:end);
mesh2_n = VHNA2.mesh{1}{2}.points(1:3);

%%

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

% Infor needed for plotting
% x1_plot_1D = [t1_mid(1): 0.001: t1_mid(end)];
% x1_plot = linspace(x1_col(1), x1_col(end), length(x1_plot_1D)).';
% y1_plot = linspace(y1_col(1), y1_col(end), length(x1_plot_1D)).'; 
% 
% x2_plot_1D = [t2_mid(1): 0.001: t2_mid(end)];
% 
% % x2_plot_1D = linspace(t2_mid(1), t2_mid(1), 1e4).';
% x2_plot = linspace(x2_col(1), x2_col(end), length(x2_plot_1D)).';
% y2_plot = linspace(y2_col(1), y2_col(end), length(x2_plot_1D)).';
% R = 20; 
%%
N_approx = 2^-6;

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
plot(mesh1_p/L1, -40*ones(length(mesh1_p)), 'o', 'DisplayName', '$\mathcal{M}_{1}^{+}$ mesh points')

hold on
plot(mesh1_n/L1, -40*ones(length(mesh1_n)), '*', 'DisplayName', '$\mathcal{M}_{1}^{-}$ mesh points')
% label('$\mathcal{M}_{1}^{-}$ mesh points')

title('Approximation of $\phi_{1}^{(2r)}$', 'fontsize',18,'interpreter','latex')
xlabel('$s/L$', 'fontsize',18,'interpreter','latex')
ylabel('$\phi_{1}^{(2r)}(s)$', 'fontsize',18,'interpreter','latex')
xlim([-0.05 1.05])
legend show

subplot(2, 1, 2)
plot(mesh2_p/L2, -40*ones(length(mesh2_p)), 'o', 'DisplayName', '$\mathcal{M}_{2}^{+}$ mesh points')
hold on
plot(mesh2_n/L2, -40*ones(length(mesh2_n)), '*', 'DisplayName', '$\mathcal{M}_{2}^{-}$ mesh points')

title('Approximation of $\phi_{2}^{(2r+1)}$', 'fontsize',18,'interpreter','latex')
xlabel('$s/L$', 'fontsize',18,'interpreter','latex')
ylabel('$\phi_{2}^{(2r+1)}(s)$', 'fontsize',18,'interpreter','latex')
xlim([-0.05, 1.05])
legend show


