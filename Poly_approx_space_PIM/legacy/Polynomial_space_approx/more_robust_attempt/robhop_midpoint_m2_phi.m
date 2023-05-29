function I = robhop_midpoint_m2_phi(x_int, y_int, x_col, y_col, FUN_m2, h, N)
% robust edditing attempt
%
% This is the midpoint function which will evaluate the integral of FUN(x)
% between a and b (where support is)
%
% Problem parameters:
% (x_int, y_int) are the coordinates of the function integration nodes (midpoints in this case)
% (x_col, y_col) are the coordinates for the collocation points
%
% OLD I THINK
% a is start of the integral
% b is the end of the integral
% FUN_m2 is the function m2(s) = k(s) - sigma1(s)m1_tilde(s)
% FUN_phi is the basis function, piecewise constants
% x is the fixed point we are evlauating the integral at
%
% Discretisation parameters:
% N is the numbe rof quadrature points we are evaluating this integral with
%
if N ~= 1
    error('Q needs to be equal to 1 for this routine to work')
end
Q = N+1;  % this gets around a blow up query!

h_new = h/Q;

% here we are taking the integration variable back to the start of the 
% interval then moving in by the new h measurement
x_nodes = [x_int - h/2 + h_new/2: h_new: x_int + h/2 - h_new/2];
y_nodes = [y_int - h/2 + h_new/2: h_new: y_int + h/2 - h_new/2];
%%% ARE THESE NODES WRONG???? I don't like it, it seems clumsey


I = 0;  % intialising integral

for j = 1:Q
    
    dist = sqrt( (x_col - x_nodes(j))^2 + (y_col - y_nodes(j))^2 );
    
    I = I + h_new*FUN_m2( dist ); %%*FUN_phi(a, b, n(j)); - REMOVED as not needed and doesn't make sense
    
%     if FUN_phi(a, b, n(j)) ==0
%         keyboard
%     else
%     end
end

end
% z = sum(w*FUN_m2(x - n)*FUN_phi(a, b, n));