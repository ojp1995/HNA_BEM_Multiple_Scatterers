function A = single_screen_mat(k, x_col, t_mid, t, C1, C2, f, h)
% In this function we will compute matrix for the single screen case.
% Problem parameters
% k is the wavenumber
% x_col is the collocation points, will be looping over this
% t_mid is the midpoints (integration nodes)
% t is the grid points
% C1 and C2 are constants for the smoothing function
% f is the function we are multiplying the hankel function in the integral
% by
%
% Discretise paramters:
% h is the step size

for j = 1: length(x_col)% looping over each collocation point
    
    
% % % %     select = (t_mid ~= x_col(j));
    A(:, j) = nosumPIM_mid_hankel(k, x_col(j), t_mid, t, C1, C2, f, h);
    
    % should this be A(j, :) so it matches with the S12_operator.m?
    
    % For each collocation point we are computing a vector that will serve
    % as a column corresponding to integating over each onterval at the
    % specific collocation point.
    
end

end