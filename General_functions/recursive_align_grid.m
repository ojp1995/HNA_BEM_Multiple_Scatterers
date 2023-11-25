function grid_quad = recursive_align_grid(grid_bf_point, grid_quad, tol)
% This is the function that will add the grid points
% WARNING: Recursive, has a get out clause.


select1 = (grid_bf_point == grid_quad);  % case where there is exact alignment already
select2 = abs(grid_bf_point - grid_quad) < tol;

if sum(select1) == 1 % this is the rare case where nothing needs doing
    disp('everything all good, no change')

elseif sum(select2) == 1  % case we need to add a point
    % should point be added before or after the select point?
    if (grid_bf_point > grid_quad(select2))
        index = find(select2);
        grid_quad = [grid_quad(1:index); grid_bf_point; grid_quad(index+1:end)];

    else
        index = find(select2);
        grid_quad = [grid_quad(1:index-1); grid_bf_point; grid_quad(index:end)];


    end    

% elseif sum(select2) > 2  % restart with smaller tolerance
%     if tol/2 < 1e-300
%         error('tolerance too small, something gone wrong')
%     end
elseif sum(select2) > 1
    % what if we just try to find the closest one?
    diff = abs(grid_bf_point - grid_quad);
    index = find(diff == min(diff));

    if grid_bf_point > grid_quad(index)
        grid_quad = [grid_quad(1:index); grid_bf_point; grid_quad(index+1:end)];

    else
        grid_quad = [grid_quad(1:index-1); grid_bf_point; grid_quad(index:end)];
%     grid_quad = recursive_align_grid(grid_bf_point, grid_quad, tol/2);
    end

end






