function [t_grid, t_mid, w, Q] = get_graded_midpoint_quad_points(L, Lgrad, ...
    h, alpha)
% In this function we will compute the graded mesh and associated weights
% for the a grid with grading towards each end a uniform mesh in the
% middle.
%
% Inputs:
%
% Outputs:
% t_grid - discretisation grid
% t_mid, midpoints for the discretisation grid
% w - length of each interval.

% computing number of points in the graded parts
% Q1 = ceil(1/( 1 - (1 - h/Lgrad)^(-alpha)) );
Q1 = ceil(alpha*Lgrad/h);

% number of quadrature points for central part
Q2 = ceil((L - 2*Lgrad)/h);

Q = 2*Q1 + Q2 ;  % number of quadrature points
t_grid = zeros(Q+1, 1);

for j = 1:Q1
    t_grid(j) = Lgrad*((j - 1)/Q1)^alpha;

    t_grid(end - j + 1) = L - t_grid(j);
end

for j = 0: Q2 
    t_grid(Q1 + j + 1) = Lgrad + j*h;
end

% computing midpoints and weights.
t_mid  = (t_grid(2:end) + t_grid(1:end-1))/2;
w = t_grid(2:end) - t_grid(1:end-1);
