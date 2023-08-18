function [t_grid, t_mid, w, h] = get_graded_midpoint_secondhalf_interval(L, ...
    Lgrad, Q2, alpha)
% second half of the interval
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
t_grid(1) = L/2;
for j = 1: Q2-1 
    t_grid(j+1) = L/2 + j*h;
end

for j = 0:Q1
    t_grid(Q1 + Q2 - j + 1) =  L - Lgrad*((j)/Q1)^alpha;
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