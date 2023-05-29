function [A, x_col, A_sing, A_smooth] = nested_int_Rob_hop_screen_mat_poly_approx_PIM(x_int, y_int, x_col, y_col, h, k, C1, C2, Q, t_disc)
% adapted to use three general integrators, rather than specific ones.
%
% ATTEMPT TO MAKE FUNCTION MORE ROBUST by using Simons style of coding.
% this function will produce the matrix for the single screen problem for
% screen, G, with start and end points of the interval a, b and N numbers
% of dof
%
%
% IS NOT TOPLITZ, JUST RUNNING SLOWLY TO START
%
%
% Problem parameters:
% x is the vector of midpoints in the x direction
% y is a vector of midpoints in the y direction
% h is the step size
% k is the wavenumber
%
%
% Discretisation paramters
% N is the number of intervals on the screen.
% Q is the number of quadrature points used
% C1 and C2 are constants needed for the smoothing function in m1_tilde
% t is the discretisation of the parameterised screen
%
% Assumption: Square system
% 


% need to define our input functions
m1 = @(t) m1_tilde( k, t, C1, C2);

m2tilde = @(t) m2(k, t, C1, C2);
% phi = @(a, b, t) bf_midpoint_rule;

for j = 1:length(x_int)  % basis function loop    
    for l = 1: length(x_col)  % collocation loop
        
        A_sing(l, j) = int_sing_samescreen_general_f(t_disc(j), t_disc(j+1), x_int(j), y_int(j), x_col(l), y_col(l), m1, @identity_function, k, (t_disc(l)+t_disc(l+1))/2,  (t_disc(l)+t_disc(l+1))/2); %singular integrator, approximating integral of m1(s)sigma1(s)f(s)
        
        A_smooth(l, j) = int_smooth_same_screen_general_f(x_int(j), y_int(j), x_col(l), y_col(l), m2tilde, @identity_function, (t_disc(l)+t_disc(l+1))/2, h(j), Q );  % composite midpoint rule, approximating integral of m2sigma2 f

    end
end

A = A_sing + A_smooth;

end