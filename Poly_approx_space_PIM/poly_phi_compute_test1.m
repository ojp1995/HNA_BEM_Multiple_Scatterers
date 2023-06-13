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

savefile = 'polycode_test1_R_20.mat';
save(savefile, 'phi1_r', 'phi2_r', 'R', 'G1', 'G2', 'C_wl', 'k', 'theta');

