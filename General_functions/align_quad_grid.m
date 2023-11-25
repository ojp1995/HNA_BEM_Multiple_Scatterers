function [grid_quad] = align_quad_grid(grid_bf, grid_quad, tol)
% In this function we will be aligning the basis function grid and the
% quadrature grid so that the integration is exactly inbetween the support
% of the basis functions.
%
% Attempting to limit as much of the recursive nature as possible by
% introducing a get out if the tol drops below a certain level.
%
% 



for j = 2:(length(grid_bf)-1)  % looping through each basis function supports
%     select1 = (grid_bf(j) == grid_quad);  % case where there is exact alignment already
%     select2 = abs(grid_bf(j) - grid_quad) < tol;

    new_grid_quad = recursive_align_grid(grid_bf(j), grid_quad, tol);

    grid_quad = new_grid_quad;


end