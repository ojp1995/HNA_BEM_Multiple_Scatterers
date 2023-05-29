function S12 = S12_operator(k, sx, sy, tx, ty, t_1D, f, h)
% in this function we will be computing the operator that will map the
% solution of a function that lives on \Gamma_{2} onto \Gamma_{1}, we will
% be integrating the following
% I(f, x) = \int_{\Gamma_{2}} H_{0}^{(1)}(k|x - y|)f(y) ds (y) for x \in
% \Gamma_{1} and for some f.
%
% Problem parameters:
% k is the wavenumber
% (sx, sy) is the coordinates for the collocation points we will be looping
% over (vector)
% (tx, ty) are the integration node midpoints 
% t_1D is the parameterised midpoints
% f is the general function in the integral
%
% Discretisation parameters:
% h is the step size (weights for integrals)

for j = 1:length(sx)
    
%     S12(:, j) = nosummidpoint_hankel_f_2D(k, sx(j), sy(j), tx, ty, t_1D, h, f);
    S12(j, :) = nosummidpoint_hankel_f_2D(k, sx(j), sy(j), tx, ty, t_1D, h, f);
    
    
end

end