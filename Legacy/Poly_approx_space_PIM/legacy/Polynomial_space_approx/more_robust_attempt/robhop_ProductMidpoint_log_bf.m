function I = robhop_ProductMidpoint_log_bf(a, b, t_col_disc, x_int, y_int, x_col, y_col, m, k, h)
%robust editting attempt
%
% We are integrating I(x) = \int_{a}^{b} m(x - t) ln|x - t| \phi(t) dt using the
% product integration rule with the midpoint method, we are approximating
% the integral I(x) as I(x) \approx \sum_{j = 1}^{N} m(x - n_{j}) \phi(n_j)
% \int_{x_{j}}^{x_{j+1}} ln(x - t) d t, for some x either in the interval
% or not.
% 
% Have added in a line for the support of the basis functions
%
% Problem parameters:
% a is the start of the interval
% b is then end point of the interval
% (x_int, y_int) are the coordinates of the function integration nodes (midpoints in this case)
% (x_col, y_col) are the coordinates for the collocation points
% m is the smooth function, m(x, t)
%
% Discretisation paraemters:
% h is the step size 
%
% Assumption:
% 1. x1 and x2 are scalar
% 2. may need to look at a different formula for when x is far away from
% [a, b].
% 3. Using midpoint integration, can only take one input node as well





dist = sqrt( (x_col - x_int)^2 + (y_col - y_int)^2 );



% this is now fixed
I = m( dist )*v1_weights_mid(a, b, k, t_col_disc); % ??? needs to be the parameterised collocation points

% old kernel idea: *v1_weights_mid( x - h/2, x + h/2, k, x);
%v1_weights_mid might have wrong inputs here
end