function z = bf_midpoint_rule(a, b, x)
% This is the basis function for the polynomial approximation where the
% basis function is a piecewise constant, taken at the middle of the
% interval. Effectively this is an indicator function
%
% Problem parameters:
% a is the begining of the support
% b is the end of the support
% x is the point we are evaluating at


if (x>=a) && (x<b)
    z = 1;
else
    z = 0;
end