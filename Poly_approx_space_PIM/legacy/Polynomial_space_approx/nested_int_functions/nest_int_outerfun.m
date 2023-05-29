% this will be the outer function that computes the scattered field for our
% iterative method. However, we will be computing each point of each matrix
% individually using a 3 different functions that integrate the three
% different integrals with a different function fed in as a variable, there
% will be nested integrals at some point (maybe not in this polynomial 
% approximation space) so these functions will become highly nested.

clear all
len = 2  % length of the screen
d = 2*pi; % distance between the screensx
G1 = [-2*pi, 2*pi, 0, 0];
G2 = [d - len*pi, len*pi, d, 0]; 

C_wl= 1/20

k = 5;  % wavenumber

theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

R_it = 10;  % R_it = 13 (I think) for len =0.5. R_it > 50, len = 2;
R_true = 20;

[err_normG1_it, err_normG2_it, aj_1_r, aj_2_r] = nested_int_it_conv_test_fixed_dof( G1, G2, k, C1, C2, theta, C_wl, R_it, R_true); 

% computing variables for G1:
[x1, y1, t1, h1, h1vector, N1, L1] = discretisation_variables(G1, C_wl, k);
% computing variables for G2:
[x2, y2, t2, h2, h2vector, N2, L2] = discretisation_variables(G2, C_wl, k);

%%
figure(1)
plot(((1:N1) - 0.5)/N1, real(aj_1_r(:, end)), 'DisplayName', '"True solution"')
hold on

figure(2)
plot(((1:N2) - 0.5)/N2, real(aj_2_r(:, end)), 'DisplayName', '"True solution"')
hold on

for r = 1: R_it
    
    figure(1)
    plot(((1:N1) - 0.5)/N1, real(aj_1_r(:, r)), 'DisplayName', sprintf('r = %g', 2*r - 2))
    
    figure(2)
    plot(((1:N2) - 0.5)/N2, real(aj_2_r(:, r)), 'DisplayName', sprintf('r = %g', 2*r-1)) 
    
end

figure(1)
legend show
title('Solution on $\Gamma_{1}$ for increasing number of iterations')
xlabel('$x/L$')
ylabel('$\phi_{1}^{r}$')
xlim([-0.1 1.1])



figure(2)
title('Solution on $\Gamma_{2}$ for increasing number of iterations')
xlabel('$x/L$')
ylabel('$\phi_{2}^{r}$')
xlim([-0.1 1.1])
legend show