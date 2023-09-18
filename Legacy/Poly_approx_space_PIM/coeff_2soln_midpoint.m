function [z, x_sample] = coeff_2soln_midpoint(aj, L, N_sample, N)
% In this function we loop over every value of x to get the solution

% Problem parameters

% Discretisation parameters


x_sample = [0: L/N_sample: L];
z = zeros(length(x_sample), 1);

for j = 1:length(x_sample)  % loopin over all the points we are sampling the solution at
    
    z(j, 1) = coeff_2_soln_midpoint_individual(aj, L, x_sample(j), N);
    
end

end