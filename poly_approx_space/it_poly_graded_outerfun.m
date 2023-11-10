% outer function for iterative polynomial solver, graded PIM and graded
% midpoint being used.

clear all

% only adding one path to centralise solvers
addpath('../General_functions/')

% introducing screens
G1_data.G = [-2*pi, 2*pi, 0, 0];

G2_data.G = [2*pi, 0, 5*pi, 3*pi]; 

Lgrad_coeff = 0.15;
alpha = 2;

C_wl= 1/20;

k = 10;  % wavenumber

theta = 0;

% constants needed for the smoothing function
C1 = 1;
C2 = pi;

% The support for the basis functions
C_wl_bf1 = 1/10;
C_wl_bf2 = 1/10;
% basis function and collocation grid
[G1_data.x_1_col, G1_data.y_1_col, G1_data.x_2_col, G1_data.y_2_col,...
    G1_data.t_bf_grid, G1_data.t_mid_col, ~, ~, G1_data.L] = discretistion_vars_graded(...
    G1_data.G, C_wl_bf1, k, Lgrad_coeff, alpha);

[G2_data.x_1_col, G2_data.y_1_col, G2_data.x_2_col, G2_data.y_2_col,...
    G2_data.t_bf_grid, G2_data.t_mid_col, ~, ~, G2_data.L] = discretistion_vars_graded(...
    G2_data.G, C_wl_bf2, k, Lgrad_coeff, alpha);

G1_data.s = [ G1_data.t_mid_col(1:end) ; flip(G1_data.L - ...
    G1_data.t_mid_col(1:end)) ];
G1_data.x_col = [ G1_data.x_1_col(1:end) ; flip(G1_data.x_2_col(1:end)) ];
G1_data.y_col = [ G1_data.y_1_col(1:end) ; flip(G1_data.y_2_col(1:end)) ];

% col_choice2 = sort(randi(length(G2_data.t_mid(:)), 20, 1));
G2_data.s = [ G2_data.t_mid_col(1:end) ; flip(G2_data.L - ...
    G2_data.t_mid_col(1:end)) ];
G2_data.x_col = [ G2_data.x_1_col(1:end) ; flip(G2_data.x_2_col(1:end)) ];
G2_data.y_col = [ G2_data.y_1_col(1:end) ; flip(G2_data.y_2_col(1:end)) ];


% quadrature nodes and other information needed
[G1_data.x_1_q, G1_data.y_1_q, G1_data.x_2_q, G1_data.y_2_q, ...
    G1_data.t_grid, G1_data.t_mid_q, G1_data.w, G1_data.N, G1_data.L] = ...
    discretistion_vars_graded(G1_data.G, C_wl, k, Lgrad_coeff, alpha);

[G2_data.x_1_q, G2_data.y_1_q, G2_data.x_2_q, G2_data.y_2_q, ...
    G2_data.t_grid, G2_data.t_mid_q, G2_data.w, G2_data.N, G2_data.L] = ...
    discretistion_vars_graded(G2_data.G, C_wl, k, Lgrad_coeff, alpha);


[S11, S12, S21, S22, u_inc1, u_inc2] = ...
    compute_matrices_for_iterative_solve(G1_data, G2_data, k, ...
    G1_data.t_bf_grid, G2_data.t_bf_grid, theta, C1, C2 );

% iterative solve
R_max = 20;

[aj_1_R, aj_2_R, phi_1_r, phi_2_r] = iterative_poly_graded_PIM_solve(...
    S11, S12, S21, S22, u_inc1, u_inc2, R_max, G1_data.t_bf_grid, ...
    G2_data.t_bf_grid, G1_data, G2_data);

%%
figure();
for r = 1:R_max
    txt1 = ['r = ', mat2str(2*r-2)];
    plot(G1_data.s/G1_data.L, phi_1_r(:, r), 'DisplayName', txt1);
    hold on

end
legend show
xlabel('$x/L_{1}$')
ylabel('$\phi_{1}^{(r)}$')
title('Iterative approximation to $\phi_{1}$ ')

figure();
for r = 1:R_max
    txt1 = ['r = ', mat2str(2*r - 1)];
    plot(G2_data.s/G2_data.L, phi_2_r(:, r), 'DisplayName', txt1);
    hold on

end
legend show
xlabel('$x/L_{2}$')
ylabel('$\phi_{2}^{(r)}$')
title('Iterative approximation to $\phi_{2}$ ')

%% solution in the domain

x_min = min([G1_data.G(1), G1_data.G(3), G2_data.G(1), G2_data.G(3)])
x_max = max([G1_data.G(1), G1_data.G(3), G2_data.G(1), G2_data.G(3)])

y_min = min([G1_data.G(2), G1_data.G(4), G2_data.G(2), G2_data.G(4)])
y_max = max([G1_data.G(2), G1_data.G(4), G2_data.G(2), G2_data.G(4)])

X1 = [-5 + x_min: 0.1 :5 + x_max];

X2 = [-5 + y_min: 0.1 :5 + y_max];

[XX, YY] = meshgrid(X1, X2);

ui = incident(k, theta, XX, YY);

% collocation nodes and phi evaluated at those nodes!
phi_1 = graded_coeff_2_solution(aj_1_R(:, end), G1_data.t_bf_grid,...
    G1_data.t_mid_q, G1_data.L);

phi_2 = graded_coeff_2_solution(aj_2_R(:, end), G2_data.t_bf_grid,...
    G2_data.t_mid_q, G2_data.L);

us = soln_in_D_2_slow(G1_data, phi_1, G2_data, phi_2, k, X1, X2);

figure();
pcolor(XX, YY, real(ui - us)); shading interp

keyboard
us_G1_slow = soln_in_D_slow(G1_data, phi_1, k, X1, X2);

us_G2_slow = soln_in_D_slow(G2_data, phi_2, k, X1, X2);


us_G1_r = soln_in_D(XX, YY, [G1_data.x_1_q ; flip(G1_data.x_2_q)],...
    [G1_data.y_1_q ; flip(G1_data.y_2_q)], k, ...
    [G1_data.w; flip(G1_data.w)], phi_1);

us_G2_r = soln_in_D(XX, YY, [G2_data.x_1_q ; flip(G2_data.x_2_q)],...
    [G2_data.y_1_q ; flip(G2_data.y_2_q)], k, ...
    [G2_data.w; flip(G2_data.w)], phi_2);

figure(); pcolor(XX, YY, real(ui - us_G2_r)); shading interp

phi_1_0 = graded_coeff_2_solution(aj_1_R(:, 1), G1_data.t_bf_grid,...
    G1_data.t_mid_q, G1_data.L);

phi_2_1 = graded_coeff_2_solution(aj_2_R(:, 1), G2_data.t_bf_grid,...
    G2_data.t_mid_q, G2_data.L);

us_G1_slow_r0 = soln_in_D_slow(G1_data, phi_1_0, k, X1, X2);

uTr0 = ui - us_G1_slow_r0;

figure(); pcolor(XX, YY, real(uTr0)); shading interp
title('R = 0')

us_G2_slow_r1 = soln_in_D_slow(G2_data, phi_2_1, k, X1, X2);

uTr1 = uTr0 - us_G2_slow_r1;
figure(); pcolor(XX, YY, real(uTr1)); shading interp
title('R = 1')
% 
% us_G1_r0 = soln_in_D(XX, YY, [G1_data.x_1_q ; flip(G1_data.x_2_q)],...
%     [G1_data.y_1_q ; flip(G1_data.y_2_q)], k, ...
%     [G1_data.w; flip(G1_data.w)], phi_1_0);
% 
% us_G2_r1 = soln_in_D(XX, YY, [G2_data.x_1_q ; flip(G2_data.x_2_q)],...
%     [G2_data.y_1_q ; flip(G2_data.y_2_q)], k, ...
%     [G2_data.w; flip(G2_data.w)], phi_2_1);

phi_1_2 = graded_coeff_2_solution(aj_1_R(:, 2), G1_data.t_bf_grid,...
    G1_data.t_mid_q, G1_data.L);

phi_2_3 = graded_coeff_2_solution(aj_2_R(:, 2), G2_data.t_bf_grid,...
    G2_data.t_mid_q, G2_data.L);

us_G1_slow_r2 = soln_in_D_slow(G1_data, phi_1_2, k, X1, X2);

uTr2 = uTr1 - us_G1_slow_r2;

figure(); pcolor(XX, YY, real(uTr2)); shading interp
title('R = 2')

us_G2_slow_r3 = soln_in_D_slow(G2_data, phi_2_3, k, X1, X2);

uTr3 = uTr2 - us_G2_slow_r3;
figure(); pcolor(XX, YY, real(uTr3)); shading interp
title('R = 3')
% 
% us_G1_r2 = soln_in_D(XX, YY, [G1_data.x_1_q ; flip(G1_data.x_2_q)],...
%     [G1_data.y_1_q ; flip(G1_data.y_2_q)], k, ...
%     [G1_data.w; flip(G1_data.w)], phi_1_2);
% 
% us_G2_r3 = soln_in_D(XX, YY, [G2_data.x_1_q ; flip(G2_data.x_2_q)],...
%     [G2_data.y_1_q ; flip(G2_data.y_2_q)], k, ...
%     [G2_data.w; flip(G2_data.w)], phi_2_3);
% 
% phi_1_4 = graded_coeff_2_solution(aj_1_R(:, 3), G1_data.t_bf_grid,...
%     G1_data.t_mid_q, G1_data.L);
% 
% phi_2_5 = graded_coeff_2_solution(aj_2_R(:, 3), G2_data.t_bf_grid,...
%     G2_data.t_mid_q, G2_data.L);
% 
% us_G1_r4 = soln_in_D(XX, YY, [G1_data.x_1_q ; flip(G1_data.x_2_q)],...
%     [G1_data.y_1_q ; flip(G1_data.y_2_q)], k, ...
%     [G1_data.w; flip(G1_data.w)], phi_1_4);
% 
% us_G2_r5 = soln_in_D(XX, YY, [G2_data.x_1_q ; flip(G2_data.x_2_q)],...
%     [G2_data.y_1_q ; flip(G2_data.y_2_q)], k, ...
%     [G2_data.w; flip(G2_data.w)], phi_2_5);

ui = incident(k, theta, XX, YY);
u = ui - (us_G1_r + us_G2_r);

figure();
pcolor(XX, YY, real(u))
hold on
% G1_x1 = linspace(G1_data.G(1), G1_data.G(3));
% G1_x2 = linspace(G1_data.G(2), G1_data.G(4));
% plot(G1_x1, G1_x2, 'LineWidth', '3')
plot([G1_data.G(1),G1_data.G(3)],[ G1_data.G(2),G1_data.G(4)], 'LineWidth', 3)
plot([G2_data.G(1),G2_data.G(3)],[ G2_data.G(2),G2_data.G(4)], 'LineWidth', 3)
shading interp

figure();
pcolor(XX, YY, real(ui - us_G1_r))
shading interp


figure();
pcolor(XX, YY, real(ui + us_G2_r))
shading interp