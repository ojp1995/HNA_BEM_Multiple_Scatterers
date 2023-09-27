function [S11, S12, S21, S22, u_inc1, u_inc2] = ...
    compute_matrices_for_iterative_solve(G1_struct, G2_struct, k, ...
    t1_bf_grid, t2_bf_grid, theta, C1, C2 )
% In this function we will be computing the 4 matrices S11, S12, S21 and
% S22 to either be used for the all in one or iterative solver, the
% corresponding right hand side vector will also be given for the initial
% incident case.
%
% The matrices will be computed with piecewise constant basis functions
% where the intervals they are supported on will be given as an input.
%
% Input parameters:
% G1 - Structure for all information needed for Gamma1
% G2 - Structure for all information needed for Gamma2
% k wavenumber
% theta, angle of incident wave  measured from downwards vertical ccw 
% C1, C2, constants for product integration method
% t1_bf_grid - Grid of support for basis functions for Gamma1
% t2_bf_grid - Grid of support for basis functions for Gamma2
%
% Assumptions:
%   1. 2D paramerterised Quadrature points given/computed for second are
%   increasing, i.e., are as in the maths from L/2 -> L, not from L -> L/2.
%
%   2. Grid of support is split into half and only half is given  


% need to either compute or unpack collocation points, Grid of support, mid
% points etc and thetransformed version of either
% unpacking structure
% G1 = G1_struct.G1;
x1_col = G1_struct.x_col;
y1_col = G1_struct.y_col;
s1 = G1_struct.s;
% unpacking quadrature points
% nq1 = G1_struct.nq;
t1_mid = G1_struct.t_mid;
t1_grid = G1_struct.t_grid;
x1_1_q = G1_struct.x_1_q;
y1_1_q = G1_struct.y_1_q;
x1_2_q = G1_struct.x_2_q;
y1_2_q = G1_struct.y_2_q;
w1 = G1_struct.w;
L1 = G1_struct.L;

% G2 = G2_struct.G2;
x2_col = G2_struct.x_col;
y2_col = G2_struct.y_col;
s2 = G2_struct.s;
% unpacking quadrature points
% nq2 = G2_struct.nq2;
t2_mid = G2_struct.t_mid;
t2_grid = G2_struct.t_grid;
x2_1_q = G2_struct.x_1_q;
y2_1_q = G2_struct.y_1_q;
x2_2_q = G2_struct.x_2_q;
y2_2_q = G2_struct.y_2_q;
w2 = G2_struct.w;
L2 = G2_struct.L;


% initialising arrays
S11 = zeros(length(x1_col), 2*length(t1_bf_grid)-2 );
S21 = zeros(length(x2_col), 2*length(t1_bf_grid) -2 );

S12 = zeros(length(x1_col), 2*length(t2_bf_grid) -2 );
S22 = zeros(length(x2_col), 2*length(t2_bf_grid) -2 );

for n = 1:length(t1_bf_grid)- 1  % basis function loop
    % NOTE, currently this is just being coded as midpoint evaluated at
    % each basis function rather than a composite midpoint, in theory it
    % shouldn't be that hard to extend it to being a composite rule but at
    % each loop of this interval a new set of quadrature points would need
    % to be computed. 

    % computation of S11
    % computing at each collocation point in one step and working though
    % each basis function

    % Finding out which quadrature points are supported by specific basis
    % function
    select1 = (t1_bf_grid(n) <= t1_mid ); 
    select2 = (t1_bf_grid(n+1) > t1_mid);
    select =  (select1 == select2);  
    clear select1 select2
    % THOUGHT!! Instead of using fnq = select, could we use w1(select) etc?
    % first half due to way it is constructed
    S11(:, n) = graded_PIM_int_hankel_f(k, s1, w1(select), ...
        t1_mid(select), 1, [t1_grid(select); t1_grid(sum(select)+1)], C1, C2);

    % computing second half
    S11(:, 2*length(t1_bf_grid)-n-1) = graded_PIM_int_hankel_f(k, L1 - s1, ...
        w1(select), t1_mid(select), 1, [t1_grid(select); t1_grid(sum(select)+1)], C1, C2);

    S21(:, n) = midpoint_hankel_f_diff_screen(k, x2_col, y2_col, x1_1_q(select),...
        y1_1_q(select), w1(select), 1);
    
    % this needs separate quadrature points, can't reuse as in S11 case.
    % Select may also need to be changed.
    xq = flip(x1_2_q); yq = flip(y1_2_q);
    S21(:, 2*length(t1_bf_grid)-n-1) = midpoint_hankel_f_diff_screen(k, ...
        x2_col, y2_col, xq(select), yq(select), w1(select), 1);

    
    % NOT tested, old, not composite rule
    % first half due to way it is constructed
%     S11(:, n) = graded_PIM_int_hankel_f(k, s1, w1(n), t1_mid_1(n), 1, ...
%         t1_grid_1(n:n+1), C1, C2);
%     
%     % computing second half
%     S11(:, length(t1_grid)+1+n) = graded_PIM_int_hankel_f(k, L1 - s1, ...
%         w1(n), t1_mid_2(n), 1, t1_grid_2(n:n+1), C1, C2);
% 
%     S21(:, n) = midpoint_hankel_f_diff_screen(k, x2_col, y2_col, x1_1_q,...
%         y1_1_q, w1(n), 1);
% 
%     S21(:, length(t1_grid)+1+n) = midpoint_hankel_f_diff_screen(k, ...
%         x2_col, y2_col, x1_2_q, y1_2_q, w1(n), 1);
    % SLOW APPROACH MAY BE EASIER Might be needed if we wanted to take a
    % composite approach but would
    % become very computationally expensive
%     for m = 1:length(x1_col)  % collocation loop
% 
%         S11(m, n) = ;  % holding for structure
% 
%     end
% 
%     % computation of S21
%     for m = 1:length(x2_col) % collocation loop for S21
% 
%         S21(m, n) = 0;  % holding for structure
% 
%     end 

end



for n = 1:length(t2_bf_grid) - 1  % basis function loop

    select1 = (t2_bf_grid(n) <= t2_mid ); 
    select2 = (t2_bf_grid(n+1) > t2_mid);
    select =  (select1 == select2);  
    clear select1 select2

        % first half
    S22(:, n) = graded_PIM_int_hankel_f(k, s2, w2(select), t2_mid(select), 1, ...
        [t2_grid(select); t2_grid(sum(select)+1)], C1, C2);

     % second half
    S22(:, 2*length(t2_bf_grid)-n - 1) = graded_PIM_int_hankel_f(k, L2 - s2, ...
        w2(select), t2_mid(select), 1, [t2_grid(select); t2_grid(sum(select)+1)], C1, C2);

%     first half
    S12(:, n) = midpoint_hankel_f_diff_screen(k, x1_col, y1_col, x2_1_q(select),...
        y2_1_q(select), w2(select), 1);
    
        % second half - As above, this is a bit of a special case, will
        % need more thought.
    xq = flip(x2_2_q); yq = flip(y2_2_q);
    S12(:, 2*length(t2_bf_grid)-n - 1) =  midpoint_hankel_f_diff_screen(k, ...
        x1_col, y1_col, xq(select), yq(select), w2(select), 1);



%     % see note on fact this is technically standard midpoint, albeit on a
%     % graded mesh, rather than a composite graded mesh.
%     % first half
%     S22(:, n) = graded_PIM_int_hankel_f(k, s2, w2(n), t2_mid_1(n), 1, ...
%         t2_grid_1(n:n+1), C1, C2);
% 
%     % second half
%     S22(:, length(t2_grid)+1+n) = graded_PIM_int_hankel_f(k, L2 - s2, ...
%         w2(n), t2_mid_2(n), 1, t2_grid_2(n:n+1), C1, C2);
% 
%     % first half
%     S12(:, n) = midpoint_hankel_f_diff_screen(k, x1_col, y1_col, x2_1_q,...
%         y2_1_q, w2(n), 1);
% 
%     % second half
%     S12(:, length(t2_grid)+1+n) =  midpoint_hankel_f_diff_screen(k, ...
%         x1_col, y1_col, x2_2_q, y2_2_q, w2(n), 1);
%    
    % SLOW - although more straight forward, here 
    % computation of S22
%     for m = 1:length(x2_col)  % collocation loop
% 
%         S22(m, n) = 1;  % holding for structure
% 
%     end
% 
%     % computation of S12
%     for m = 1:length(x1_col) % collocation loop for S12
% 
%         S12(m, n) = 0;  % holding for structure
% 
%     end 

end

% computing incident for each

u_inc1 = incident(k, theta, x1_col, y1_col);
u_inc2 = incident(k, theta, x2_col, y2_col);





