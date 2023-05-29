function z = coeff_2_soln_midpoint_individual(aj, L, x, N)
% Will take the coefficients and compute the solution

% Problem parameters

% Discretisation parameters


h = L/N; %step size
t = [0: h: L]; % discretisation of screen for basis functions
% looping over the interval
z = 0;
for j = 1: N

    z = z + aj(j)*bf_midpoint_rule(t(j), t(j+1), x);
        
end