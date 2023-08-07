function phi_j_r = get_phi_j_r(coeffs, Gj, Lj, k, d, h_ell, y1_ell, ...
    y2_ell, nj, t1d, xj, yj, phi_jn1_rn1, alpha)
% In this function we will compute the solution phi_j^r from the
% coefficients and leading order behaviour at the points we want it
% evaluated at.
%
% Inputs
% coeffs, these are the coefficients
% LoB, the leading order behaviour
% phi_j-1_^r-1, the previous phi 
%
% t1d, parameterised point to evaluate at
% Gj, Gamma_j, we are evlauting the solution on
% Lj, length of Gamma_j
% k, the wavenumber
% d, direction of incident wave
% alpha factor for uinc LoB

LoB = midpoint_dphikdn_f_diff_screen(k, xj, yj, h_ell, y1_ell, y2_ell, ...
    phi_jn1_rn1.', nj);

phi_j_r = coeffs.eval(t1d, 1) + 2*alpha*duidn(Gj, Lj, k, d, t1d).' + LoB;

% LoB sign is correct, it is a minus until the derivative comes into play
% as the d/dx (H_0^(1)(x)) = -H_1^(1)(x) giving the plus overall.