function [x_1, y_1, x_2, y_2, t_grid, t_mid, w, N, L] = ...
    discretistion_vars_graded(G, C_wl, k, Lgrad_coeff, alpha)
% In this function we will compute the graded grid for the half
% interval and the midpoints in the parameterised and in 2D for the half
% interval as well as a vector of the weights for the half interval, the
% number of discretisation points and the length of the screen.

% first computing the length of the screen
L = sqrt( (G(3) - G(1))^2 + (G(4) - G(2))^2 ); 

Lgrad = L*Lgrad_coeff;

% number of dof
N = ceil(k*L./(C_wl*2*pi));

% nodes and grid for half interval
[t_grid, t_mid, w, h] = get_graded_midpoint_half_interval(L, ...
    Lgrad, N, alpha);

N = length(t_mid);
% computing the midpoints coordinates for \Gamma, on first half of interval
x_1 = G(1) + t_mid*(G(3) - G(1))/L;
y_1 = G(2) + t_mid*(G(4) - G(2))/L;

% computing the midpoints coordinates for \Gamma, on second half of interval
x_2 = G(1) + (L - t_mid)*(G(3) - G(1))/L;
y_2 = G(2) + (L - t_mid)*(G(4) - G(2))/L;

end