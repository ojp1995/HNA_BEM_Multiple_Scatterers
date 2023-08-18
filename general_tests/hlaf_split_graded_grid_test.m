% test to make sure that half grid method is doing what we expect it to
clear all

addpath('../General_functions/')

L = 2;
Lgrad = 0.15*L;
Q2 = 5;
alpha = 2;

% left hand side grid:
[t_gridLHS, t_midLHS, wLHS, hLHS] =...
            get_graded_midpoint_half_interval(L, Lgrad, Q2, alpha);
t_gridLHS
% right hand side grid:
[t_gridRHS, t_midRHS, wRHS, hRHS] = get_graded_midpoint_secondhalf_interval(L, ...
    Lgrad, Q2, alpha);