function [f_ell, f_ell_uinc, f_ell_beam, ddn_S_ell_j_phij,S_ell_j_phi_j,...
    LoB_uinc, LoB_beam, LoB_kernel_uinc, LoB_kernel_beam] ...
    = beam_compute_RHS_vec_given_coll_vec(Gamma_ell, L_ell, k, d, ...
    theta, norm_j, s_ell, x1s_ell, x2s_ell, C1, C2, y1nq_j, y2nq_j, ...
    h_j_outer, phi_j_outer, y1nq_ell_PIM, y2nq_ell_PIM, nq_ell_PIM, ...
    tq_ell_PIM, h_ell_outer_PIM, y1nq_j_inner, y2nq_j_inner, h_j_inner, ...
    phi_j_inner)
% In this function ...
%
% Problem parameters:
% Gamma_ell, set up of the screen Gamma_ell
% L_ell, the length of the screen Gamma_ell
% k the wave number
% d the unit direction of the incident wave
% theta the angle of the incident wave measured from the downwards
% vertical anti-clockwise 
% norm_j, unit normal of Gamma_j
% s_ell, parameterised collocation points in 1D
% (x1s_ell, x2s_ell), Collocation points in 2D
% C1, C2, constants needed for PIM method
%
% Discretisation parameters:
% 1. Integration variables for approximating outer integral over Gamma_j:
%     (y1nq_j, y2nq_j), integration nodes for this integral
%     h_j_outer, integration weights for this integral
%     phi_j_outer, vector phi_j evlauated at the integration nodes
%
% 2. Integration variables for approximating the integral over Gamma_ell
% using the product integration method
%     (y1nq_ell_PIM, y2nq_ell_PIM), integration nodes for this integral
%     nq_ell_PIM, parameterised 1D integration nodes for this integral
%     tq_ell_PIM, discretisation on Gamma_ell for integration
%     h_ell_outer_PIM, intergation weight
%
% 3. Integration variables for approximating the integral Psi_ell_r:
%     (y1nq_j_inner, y2nq_j_inner), intergation nodes in 2D
%     h_j_inner, integration weights
%     phi_j_inner, phi_j evaluated at the integration nodes
%
% Outputs:
% f_ell, right hand side evaluated at the collocation points, s.

% First we will compute the leading order behaviour where the collocation
% points we are evaluating the integral are the quadrature nodes for the
% outer integral.

ddn_S_ell_j_phij = midpoint_dphikdn_f_diff_screen(k, y1nq_ell_PIM, ...
    y2nq_ell_PIM, h_j_inner, y1nq_j_inner, y2nq_j_inner, ...
    phi_j_inner.', norm_j).';

LoB_kernel_uinc = 2*duidn(Gamma_ell, L_ell, k, d, nq_ell_PIM);
LoB_kernel_beam = ddn_S_ell_j_phij;

% now we will compute the whole vector
% beam from Gamma_j incident onto Gamma_ell
S_ell_j_phi_j = midpoint_hankel_f_diff_screen(k, x1s_ell, x2s_ell, ...
y1nq_j, y2nq_j, h_j_outer, phi_j_outer.');

LoB_uinc =  PIM_int_hankel_f(k, s_ell, h_ell_outer_PIM, ...
nq_ell_PIM, LoB_kernel_uinc, tq_ell_PIM, C1, C2);

LoB_beam = PIM_int_hankel_f(k, s_ell, h_ell_outer_PIM, ...
nq_ell_PIM, LoB_kernel_beam, tq_ell_PIM, C1, C2);

f_ell_uinc = incident(k, theta, x1s_ell, x2s_ell) - LoB_uinc;

f_ell_beam = - S_ell_j_phi_j + LoB_beam;

f_ell = f_ell_uinc + f_ell_beam;

end

