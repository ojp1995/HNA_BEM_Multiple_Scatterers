function S21 = S21_op_coeff_robhop(x_int, y_int, x_col, y_col, k, FUN, n, h)
% Rewirtten so it evaluates the integrals instead of just the matrix
%
% REWRITTEN as robust hopefully.
% This function computes S21\phi_{1}^{(r)} = \int_{\G1}\phik(x,
% y)\phi_{1}^{(r)}(y) ds(y), x \in \G2
% Idea here is the beam originating from G1 incident on G2.
%
% Can produce S12 op if coordinates are flipped at input.
%
% Problem parameters
% (x_int, y_int) is the coordinates of G1, our integration nodes
% (x_col, y_col) are the coordinates of G2, our collocvation points
% k is the wavenumber
% FUN is the function we are evlauating, in this case the coefficients
% n is the point we are evaluating the integral at.
% Discretisation parameters:
% h are the integration weights
S21 = zeros(length(x_col), 1);
for j = 1: length(x_int)  %looping through the integration intervals on \Gamma_1
    
    for l = 1:length(x_col)  % looping through the collocation points on \Gamma_2
        S21(l, 1) = S21(l, 1) + int_diff_screen_phik_general_f(x_int(j), y_int(j), x_col(l), y_col(l), k, FUN(n(j)), h);         
    end
    
end