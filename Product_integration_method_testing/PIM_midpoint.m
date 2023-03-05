function z = PIM_midpoint(sigma1, m1, sigma2, m2, w, n)
% The aim of this function is to compute an integral which has a
% singularity in it and to get around this using the Product Integration
% Method, PIM. Using a midpoint quadrature method. The function we are
% integrating is split as k(x) = sigma1(x)m1(x) + sigma2(x)m2(x).
%
% Assumptions:
% 1. We are assumming that we are seperating the integral into two terms.
% 2. Only works for PIM as weights are equal to h.
%
% Problem parameters:
% sigma1, simple but potentially singular part of the function
% m1, smooth part 
% sigma2, simple but potentially singular part of the function
% m2, smooth part
%
% Discretisation parameters:
% w, the weights, in this case h, the step size
% n, the nodes, the midpoints of the intervals

%%%% WEIGHTS NEED TO BE SPECIFIC for each

z = midpoint_vec(m1, w, n) +  midpoint_vec(m2, w, n);

end