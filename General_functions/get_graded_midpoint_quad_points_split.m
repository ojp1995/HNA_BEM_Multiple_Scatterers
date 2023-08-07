function [t_grid, t_mid, w, Q, t_grid_end, t_mid_end, w_end] = ...
    get_graded_midpoint_quad_points_split(L, Lgrad, ...
    h, alpha)
% In this function we will compute the graded mesh and associated weights
% for the a grid with grading towards each end a uniform mesh in the
% middle.
%
% Split so the second part grading has the integral is graded towards 0 as
% well.
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
t_grid = zeros(Q1 + Q2 +1, 1);

for j = 1:Q1
    t_grid(j) = Lgrad*((j - 1)/Q1)^alpha;
    
    t_grid_end(j) = t_grid(j);
%     t_grid(end-Q1 + j) = t_grid(j);
end

for j = 0: Q2 
    t_grid(Q1 + j + 1) = Lgrad + j*h;
end

% computing the midpoints
t_mid = zeros(Q1 + Q2, 1);
for j = 1:Q1
    t_mid(j) = (t_grid(j+1) + t_grid(j))/2;
    t_mid_end(j) = t_mid(j);
%     t_mid(Q1 + Q2 + j) = (t_grid(j+1) + t_grid(j))/2;
end

for j = Q1 + 1: Q1 + Q2
    t_mid(j) = (t_grid(j+1) + t_grid(j))/2;
end
% t_mid  = (t_grid(2:end) + t_grid(1:end-1))/2;

% computing the weights
w = zeros(Q1 + Q2, 1);
for j = 1:Q1
    w(j) = t_grid(j+1) - t_grid(j);
    w_end(j) = w(j);
%     w(Q1 + Q2 + j) = w(j);
end

for j = Q1 + 1: Q1 + Q2
    w(j) = t_grid(j+1) - t_grid(j);

end
