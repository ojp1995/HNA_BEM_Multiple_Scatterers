% graded midpoint testing

clear all
addpath('../General_functions/')

L = 2;
Lgrad = 0.15*L;
h = 0.2;
alpha = 2;

[t_grid, t_mid, w, Q] = get_graded_midpoint_quad_points(L, Lgrad, h, alpha);

figure()
plot(t_grid, zeros(length(t_grid), 1), '*', 'DisplayName', 'Grid points')
hold on
plot(t_mid, zeros(length(t_mid), 1), 'o', 'DisplayName', 'Midpoints')
plot(t_grid, zeros(length(t_grid), 1), 'k-')
xlim([0 L])
legend show