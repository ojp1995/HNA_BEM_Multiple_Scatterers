function [aj_1_r, aj_2_r, A1, A2, S12, S21] = iterative_2screen_solver(k, G1, G2, f, theta, C_wl, R_approx, C1, C2)
% In this function we will be computing the solution on the boundary for 
% using the iterative method for R refelctions.
%
% The matrix A1 is the single screen problem matrix for \Gamma_{1} defined 
% as
%
% The matrix A2 is the single screen problem matrix for \Gamma_{2} defined 
% as
%
% The matrix S12 is a matrix mapping the solution on \Gamma_{2} to
% \Gamma_{1} defined as
%
% The matrix S21 is a matrix mapping the solution on \Gamma_{1} to
% \Gamma_{2} defined as
%
% Problem parameters:

[x1, y1, t1, t1_mid, h1, hvector1, N1, L1] = discretisation_variables(G1, C_wl, k);
[x2, y2, t2, t2_mid, h2, hvector2, N2, L2] = discretisation_variables(G2, C_wl, k);

% when on the same screen, collocation points and integration nodes are the
% same points.
x1_col = t1_mid;
x2_col = t2_mid;
% single screen solves

A1 = single_screen_mat(k, x1_col, t1_mid, t1, C1, C2, f, h1);

A2 = single_screen_mat(k, x2_col, t2_mid, t2, C1, C2, f, h2);

% S12 and S21 operators
S12 = S12_operator(k, x1, y1, x2, y2, t2_mid, f, h2);

S21 = S12_operator(k, x2, y2, x1, y1, t1_mid, f, h1);


% solving for solution on the boundary
% step 0
u_i_1 = robhop_PW_incident(k, theta, x1, y1).';
aj_1_r(:, 1) = A1\u_i_1;

% figure()
% plot(t1_mid, real(aj_1_r(:, 1) ))
% step 1
% error('For matrix to work should be S12 not S21 in line below.')
u_i_2 = robhop_PW_incident(k, theta, x2, y2).' - S21*aj_1_r(:, 1);  
aj_2_r(:, 1) = A2\u_i_2;

% figure()
% plot(t2_mid, real(aj_2_r(:, 1) ))
% iterative solve
for r = 2:R_approx
    
    % solve on \Gamma_{1}
    u_i_1 = robhop_PW_incident(k, theta, x1, y1).' - S12*aj_2_r(:, r-1);
    aj_1_r(:, r) = A1\u_i_1;
    
    % solve on \Gamma_{2}
    u_i_2 = robhop_PW_incident(k, theta, x2, y2).' - S21*aj_1_r(:, r);  
    aj_2_r(:, r) = A2\u_i_2;
    
    
end 

end
