function I = int_same_screnn_PIM_gen_f(a, b, x_int, y_int, x_col, y_col, m1, m2, FUN, k, h, n, t_col_disc)
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
% m1 is the smooth function, m1(x, t), with corresponding singular function
% \sigma_{1}(x, t).
% m2 is the smooth function m2(s) = k(s) - sigma1(s)m1_tilde(s)
%
% Inputs for our general function
% the point we want to evaluate it at?

% Discretisation parameters:
% w - the weights, computed in v1_weights_mid
% n - the nodes
% N is the number of intercals we are splitting this up into
% h is the step size
% Asuumptions:
% 1. only taking a scalr input for x(s).
% 2. Uses product integation method
% 3. Integration variables on the same screen as collocation point 
dist = sqrt( (x_col - x_int)^2 + (y_col - y_int)^2 );

I_sing = m1(dist)*v1_weights_mid(a, b, k, t_col_disc)*FUN(n);

if N ~= 1
    error('Q needs to be equal to 1 for this routine to work')
end
Q = N+1;  % this gets around a blow up query!

h_new = h/Q;

% here we are taking the integration variable back to the start of the 
% interval then moving in by the new h measurement
x_nodes = [x_int - h/2 + h_new/2: h_new: x_int + h/2 - h_new/2];
y_nodes = [y_int - h/2 + h_new/2: h_new: y_int + h/2 - h_new/2];
n_nodes = [n - h/2 + h_new/2: h_new: n + h/2 - h_new/2];

I_smooth = 0;  % intialising integral

for j = 1:Q
    
    dist = sqrt( (x_col - x_nodes(j))^2 + (y_col - y_nodes(j))^2 );
    
    I_smooth = I_smooth + h_new*m2( dist )*FUN(n_nodes(j)); %%*FUN_phi(a, b, n(j)); - REMOVED as not needed and doesn't make sense
    

end

I = I_sing + I_smooth;