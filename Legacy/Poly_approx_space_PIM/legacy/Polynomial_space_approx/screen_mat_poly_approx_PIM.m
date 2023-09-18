function [A, x_col, A_sing, A_smooth] = screen_mat_poly_approx_PIM(N, a, b, k, C1, C2, Q)
% this function will produce the matrix for the single screen problem for
% screen, G, with start and end points of the interval a, b and N numbers
% of dof
%
%
% IS NOT TOPLITZ, JUST RUNNING SLOWLY TO START
%
%
% Problem parameters:
% a is the start of the interval
% b is the end of the interval
% k is the wavenumber
%
%
% Discretisation paramters
% N is the number of intervals on the screen.
% Q is the number of quadrature points used
% C1 and C2 are constants needed for the smoothing function in m1_tilde
%
% Assumption: Square system
% 

% step size
h = (b - a)/N;

% collocation points, in this case we are choosing the midpoints
x_col = [a + h/2: h: b - h/2];

% need to define our input functions
m1 = @(x) m1_tilde( k, x, C1, C2);

m2tilde = @(x) m2(k, x, C1, C2);
% phi = @(a, b, t) bf_midpoint_rule;


for j = 1:N  % basis function loop    
    for l = 1: length(x_col)  % collocation loop
        
        
        A_sing(l, j) =  ProductMidpoint_log_bf(a + (j-1)*h, a + j*h, k, x_col(l), m1, @bf_midpoint_rule, Q); % singular part of the intergral
        
        A_smooth(l, j) =  midpoint_m2_phi(a + (j-1)*h, a + j*h, x_col(l), m2tilde, @bf_midpoint_rule, Q); % non-singluar part of the integral m2sigma2phi

    end

end

A = A_sing + A_smooth;
end
