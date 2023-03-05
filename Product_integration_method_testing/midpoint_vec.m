function z = midpoint_vec(FUN, w, n)
% This is the midpoint function which will evaluate the integral of FUN(x)
% between a and b.
%
% Problem parameters:
% FUN is the function we are integrating
%
% Discretisation parameters:
% w is the vector of weights
% n is the nodes, midpoint of the interval.

z = sum(w.*FUN(n));