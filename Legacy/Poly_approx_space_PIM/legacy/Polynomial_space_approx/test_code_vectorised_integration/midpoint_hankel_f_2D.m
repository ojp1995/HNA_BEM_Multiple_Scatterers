function z = midpoint_hankel_f_2D(k, sx, sy, tx, ty, t_1D, h, f)
% In this function we will be computing the midpoint rule for the following
% I(f, x) = \int_{\Gamma_{1}} H_{0}^{(1)}(k |x - y|) f(y) ds(y) for x \in
% \Gamma_{2}.
%
% Problem parameters:
% k is the wavenumber
% (sx, sy) is the position x we are evaluating the integral at (scalar)
% (tx, ty) are the integration nodes (each vectors)
% t_1D is the paramterised integration nodes (vector)
% h is the width interval/weights
% f is the function the general function in the integral that depends on
% t_1D

z = sum(1i*h*besselh(0, k*sqrt( (sx - tx).^2 + (sy - ty).^2 )).*f(t_1D)/4);