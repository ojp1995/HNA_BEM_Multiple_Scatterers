% test function for computing the vector S21phi1_0
clear all

% addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/Poly_approx_space_PIM/legacy/Polynomial_space_approx/Simon_code')
addpath('/Users/Oliver/Dropbox/Mac (2)/Documents/Github/HNA_BEM_Multiple_Scatterers/Poly_approx_space_PIM/legacy/Polynomial_space_approx/more_robust_attempt')
% introducing the screens



% switch and things get interesting
G1 = [-2*pi, 2*pi, 0, 0];

G2 = [2*pi, 0, 5*pi, 3*pi]; 

C_wl= 2^-6;

k = 10;  % wavenumber

theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

% number of iterations for approximation and truth

R_it = 20;  % R_it = 13 (I think) for len =0.5. R_it > 50, len = 2;
R_true = 20;

[err_normG1_it, err_normG2_it, aj_1_r, aj_2_r, aj1_1step, aj2_1step]...
    = it_conv_test_fixed_dof( G1, G2, k, C1, C2, theta, C_wl, R_it, R_true); 

% comouting solution

% computing variables for G1:
[x1, y1, t1, h1, h1vector, N1, L1] = discretisation_variables(G1, C_wl, k);
% computing variables for G2:
[x2, y2, t2, h2, h2vector, N2, L2] = discretisation_variables(G2, C_wl, k);

R = R_true;

phi1_r = zeros(R, length(x1));
phi2_r = zeros(R, length(x2));
for r = 1:R
    
    phi1_r(r, :) = coeff_2soln_midpoint(aj_1_r(:, r), L1, N1-1, N1);
    phi2_r(r, :) = coeff_2soln_midpoint(aj_2_r(:, r), L2, N2-1, N2);
    
end
%%
h1 = h1vector(1);
t1_mid = [t1(1) + h1/2:h1 : t1(end) - h1/2];
h2 = h2vector(1);
t2_mid = [t2(1) + h2/2:h2 : t2(end) - h2/2];
figure()


for r = 1:4
    subplot(2, 1, 1)
    plot(t1_mid/L1, real(phi1_r(r, :)), 'DisplayName', strcat('r = ', num2str(2*r - 2)))
    hold on
    
    subplot(2, 1, 2)
    plot(t2_mid/L2, real(phi2_r(r, :)), 'DisplayName', strcat('r = ', num2str(2*r - 1)))
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

% savefile = 'polycode_test1_R_20.mat';
% save(savefile, 'phi1_r', 'phi2_r', 'R', 'G1', 'G2', 'C_wl', 'k', 'theta');

