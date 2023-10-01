function [S11, u_inc] = compute_matrix_single_screen_solve_graded(...
    G_struct, k, t1_bf_grid, theta, C1, C2)
% In this function we are computing the single screen solve matrix only

% unpacking structure
% G1 = G1_struct.G1;
x1_col = G_struct.x_col;
y1_col = G_struct.y_col;
s1 = G_struct.s;
% unpacking quadrature points
% nq1 = G1_struct.nq;
t1_mid = G_struct.t_mid;
t1_grid = G_struct.t_grid;
% x1_1_q = G_struct.x_1_q;
% y1_1_q = G_struct.y_1_q;
% x1_2_q = G_struct.x_2_q;
% y1_2_q = G_struct.y_2_q;
w1 = G_struct.w;
L1 = G_struct.L;

S11 = zeros(length(x1_col), 2*length(t1_bf_grid)-2 );


for n = 1:length(t1_bf_grid)- 1  % basis function loop

    % Finding out which quadrature points are supported by specific basis
    % function
    select1 = (t1_bf_grid(n) <= t1_mid ); 
    select2 = (t1_bf_grid(n+1) > t1_mid);
    select =  (select1 == select2); 
    grid_select = find(select); 
    ii = max(grid_select); 
    grid_select(length(grid_select)+1) = ii+1;
    clear select1 select2
    % THOUGHT!! Instead of using fnq = select, could we use w1(select) etc?
    % first half due to way it is constructed
    
    S11(:, n) = graded_PIM_int_hankel_f(k, s1, w1(select), ...
        t1_mid(select), 1, t1_grid(grid_select), C1, C2);
    
    % computing second half
    S11(:, 2*length(t1_bf_grid)-n-1) = graded_PIM_int_hankel_f(k, L1 - s1, ...
        w1(select), t1_mid(select), 1, t1_grid(grid_select), C1, C2);

end

u_inc = incident(k, theta, x1_col, y1_col);