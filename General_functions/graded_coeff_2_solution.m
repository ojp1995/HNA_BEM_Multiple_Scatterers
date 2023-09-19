function phi = graded_coeff_2_solution(aj, t_bf_grid, x, L)
% Given the coefficients and the half grid of support for basis functions,
% find the value of phi at the points x, where x is in the range [0, L/2].
% This allows us to take advantage of knowing we are near 0.
phi1 = zeros(length(x), 1);
phi2 = zeros(length(x), 1);
phi = zeros(2*length(x), 1);
for j = 1:length(t_bf_grid) - 1
    
    % first see which interval the points x are in
    select1 = (t_bf_grid(j) <= x ); 
    select2 = (t_bf_grid(j+1) > x);
    select =  (select1 == select2);  
    clear select1 select2

    phi1(select) = aj(j);

    phi2(select) = aj(length(aj) - j + 1);

end

phi = [phi1(:) ; flip(phi2(:))];