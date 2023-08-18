function [t_grid, t_mid, w, h] = get_graded_midpoint_half_interval(L, ...
    Lgrad, Q2, alpha)
% In this function we will be computing the quadrature points, grid and
% weights for the interval [0, L/2], with grading on the part of the 
% interval [0, Lgrad] and uniform from [Lgrad, L/2]
%
% Inputs:
% L, length of the screen (may change to length of interval we are interested in)
% Lgrad, end point of graded quadrature points
% h, uniform section step size
% alpha grading parameter

h = (L/2 - Lgrad)/Q2;

Q1 =ceil(alpha*Lgrad/h);



Q = Q1 + Q2;
% computing the grid
t_grid = zeros(Q1 + Q2 + 1, 1);
for j = 1:Q1
    t_grid(j) =  Lgrad*((j - 1)/Q1)^alpha;
end

for j = 0: Q2 
    t_grid(Q1 + j + 1) = Lgrad + j*h;
end
% computing the midpoints
t_mid = zeros(Q1 + Q2, 1);
for j = 1:Q1
    t_mid(j) = (t_grid(j+1) + t_grid(j))/2;
end

for j = Q1 + 1: Q1 + Q2
    t_mid(j) = (t_grid(j+1) + t_grid(j))/2;
end
% computing the weights
w = zeros(Q1 + Q2, 1);
for j = 1:Q1
    w(j) = t_grid(j+1) - t_grid(j);
end

for j = Q1 + 1: Q1 + Q2
    w(j) = t_grid(j+1) - t_grid(j);

end