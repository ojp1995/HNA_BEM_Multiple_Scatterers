function I = int_sing_samescreen_general_f(a, b, x_int, y_int, x_col, y_col, m, FUN, k, n, t_col_disc)
% In this funciton we will compute the following integral
% I(x(s), f) = \int_{G1} H_{0}^{(1)}(k |x(s) - y(t)|) f(y(t)) ds(y(t)) for
% x(s) \in \G1 and for some function f.
% 
% As the integration variable and point x(s) are both in the same screen we
% will be using a product integration method to get around the singularity.
%
% Problem parameters:
% a is the start of the interval
% b is then end point of the interval
% (x_int, y_int) are the coordinates of the function integration nodes (midpoints in this case)
% (x_col, y_col) are the coordinates for the collocation points
% m is the smooth function, m(x, t)
%
% Inputs for our general function
% the point we want to evaluate it at?

% Discretisation parameters:
% w - the weights, computed in v1_weights_mid
% n - the nodes
% 
% Asuumptions:
% 1. only taking a scalr input for x(s).
% 2. Uses product integation method
% 3. Integration variables on the same screen as collocation point 

dist = sqrt( (x_col - x_int)^2 + (y_col - y_int)^2 );

I = m(dist)*v1_weights_mid(a, b, k, t_col_disc)*FUN(n);