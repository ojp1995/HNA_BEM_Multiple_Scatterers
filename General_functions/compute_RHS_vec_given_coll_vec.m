function f_ell = compute_RHS_vec_given_coll_vec(Gamma_ell, L_ell, k, ...
    d, s, x1s, x2s, nq, y1nq, y2nq, h_outer, phi_j_outer, tq, C1, C2, ...
    h_inner, y1nq_inner, y2nq_inner, phi_j_inner, nj , theta)
% In this function ...
%
% Inputs:
% Gamma_ell, is the coordinates of the Gamma_ell
% L_ell is the length of Gamma_ell
% k, the wavenumber
% s, 1D parameterised collocation points
% d is the directional unit vector of the incident wave
% (x1s, x2s), collocation points evaluated in 2D
% nq, paramterised integration nodes, outer integral
% (y1nq, y2nq), integration nodes in 2D
% h_outer, integration weight for outer integral
% phi_j_outer, phi_j evaluated at outer quadrature nodes.
% tq, mesh for integration of outer integral (PIM)
% C1, C2, constants needed for smoothing 
% 
% (y1nq, y2nq), integration nodes in 2D for outer integral
% h_inner, integration weights for inner integral
% (y1nq_inner, y2nq_inner), integration nodes for inner integral
% phi_j_inner, phi_j elvaluated at inner quadrature nodes, vector.
% nj, unit normal for Gamma_j.
% theta, the angle of the incident wave measured from the downwards
% vertical anti-clockwise 

% Outputs:
% f_ell, right hand side evaluated at the collocation points, s.

% First we will compute the leading order behaviour where the collocation
% points we are evaluating the integral are the quadrature nodes for the
% outer integral.
Psi_ell_r = 2*duidn(Gamma_ell, L_ell, k, d, nq) ...
    - 2*midpoint_dphikdn_f_diff_screen(k, y1nq, y2nq, h_inner, ...
    y1nq_inner, y2nq_inner, phi_j_inner, nj);

% now we will compute the whole vector
f_ell = incident(k, theta, x1s, x2s) ...
    - midpoint_hankel_f_diff_screen(k, x1s, x2s, y1nq, y2nq, h_outer, ...
    phi_j_outer) - PIM_int_hankel_f(k, s, h_outer, nq, Psi_ell_r, tq, ...
    C1, C2); 