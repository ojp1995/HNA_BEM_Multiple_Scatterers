function I = int_smooth_same_screen_general_f(x_int, y_int, x_col, y_col, FUN_m2, FUN, n, h, N )
% In this function we are hoping to integrate phik for a general function
% f.
%
% Problem parameters:
% (x_int, y_int) are the coordinates of the function integration nodes (midpoints in this case)
% (x_col, y_col) are the coordinates for the collocation points
% FUN_m2 is the function m2(s) = k(s) - sigma1(s)m1_tilde(s)
% FUN is the general function we are approximating
% n the point we are evaluating the function at
%
% Discretisation parameters:
% N is the number of intercals we are splitting this up into
% h is the step size
%
% Asuumptions
% see error function below

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

I = 0;  % intialising integral

for j = 1:Q
    
    dist = sqrt( (x_col - x_nodes(j))^2 + (y_col - y_nodes(j))^2 );
    
    I = I + h_new*FUN_m2( dist )*FUN(n_nodes(j)); %%*FUN_phi(a, b, n(j)); - REMOVED as not needed and doesn't make sense
    

end

end