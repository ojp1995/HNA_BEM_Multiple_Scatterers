function [A, x_col, A_sing, A_smooth] = Rob_hop_screen_mat_poly_approx_PIM(x_int, y_int, x_col, y_col, h, k, C1, C2, Q, t_disc)
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
        
        % this might be wrong!, the singular part of the matrix may
        % actually be needed at all points, need to figure that out!
            % here (x(l), y(l)) are the collocation points and (x(j), y(j))
            % are the points we are integrating with respect to?
            
            % A_sing matrix is currently wrong, mistake in computing
            % weights
            A_sing(l, j) = robhop_ProductMidpoint_log_bf(t_disc(j), t_disc(j+1), (t_disc(l)+t_disc(l+1))/2 , x_int(j), y_int(j), x_col(l), y_col(l), m1, k, h(j));  % singular part of the intergral
            
            A_smooth(l, j) =  robhop_midpoint_m2_phi(x_int(j), y_int(j), x_col(l), y_col(l), m2tilde , h(j), Q); % non-singluar part of the integral m2sigma2phi
        
    end

end

A = A_sing + A_smooth;
end
