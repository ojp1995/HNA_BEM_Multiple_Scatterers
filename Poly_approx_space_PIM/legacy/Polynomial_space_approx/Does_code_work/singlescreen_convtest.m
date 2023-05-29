% single screen convergence

% adding paths and loading data
clear all
addpath('/Users/ojp18/OneDrive - University of Reading/PhD/MATLAB/Polynomial_space_approx/') % code to be tested


G1 = [-2*pi, 2*pi, 0, 0];
L1 = sqrt( (G1(3) - G1(1))^2 + (G1(4) - G1(2))^2 );  % length of G1
G2 = [4*pi, 0, 7*pi, 3*pi];
L2 = sqrt( (G2(3) - G2(1))^2 + (G2(4) - G2(2))^2 );  % length of G2

k = 5;  % wavenumber

theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

% I think the set up has now been introduced so we can start to compute the
% matrices for \Gamma_{1} and \Gamma_{2}

Q = 1;  % quadrature discretisation

C_wl = [1/5, 1/10, 1/20, 1/40, 1/80, 1/160, 1/320]; % h <= C_wl \lambda, 1/C_wl intervals per wl.

% computing "true" solution
C_wl_true = 1/640;
N1 = ceil(k*L1./(C_wl_true*2*pi));
[A1, x_col1_true, A1_sing, A1_smooth] = screen_mat_poly_approx_PIM(N1, 0, L1, k, C1, C2, Q);
u_i_1_r = PW_incident(k, theta, G1, x_col1_true);
aj_true = A1\u_i_1_r.';

figure()
hold on
for j = 1:length(C_wl)
    r = 1;
    % number of intervals on wach screen:
     N1 = ceil(k*L1./(C_wl(j)*2*pi)) % number of itervals on G1
     
     [A1, x_col1, A1_sing, A1_smooth] = screen_mat_poly_approx_PIM(N1, 0, L1, k, C1, C2, Q);
  
    % solutions for \Gamma_{1}
    u_i_1_r = PW_incident(k, theta, G1, x_col1);
    aj_1_r = A1\u_i_1_r.';
    
    
    plot(x_col1/L1, aj_1_r, 'DisplayName', sprintf('N = %g', 1/C_wl(j)))
    
    % computiong error
    % first mapping true solution to coarse grid
    keyboard
    phi1_true = interp1(x_col1_true, aj_true, x_col1);
    
    err_normG1_it_interpl1(j) = norm(phi1_true - aj_1_r.', 1)/norm(phi1_true, 1);
    
end

legend show
title('Single screen solution, PIM method')
xlabel('x/L')
ylabel('\phi_{1, N}')

% estimated order of convergence
for j = 1:(length(C_wl) - 1)
    
    EOC(j)= log2(err_normG1_it_interpl1(j)/err_normG1_it_interpl1(j+1));
    
end


% plotting the error


