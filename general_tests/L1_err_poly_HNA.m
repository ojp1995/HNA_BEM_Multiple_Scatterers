function [err_no_norm, err_norm_poly_true, err_norm_HNA_true, phi_poly, phi_HNA] = L1_err_poly_HNA(tq, xq, yq, w, G, L, k, d, n,...
    poly_coeffs, HNA_coeffs,  h_in_int, y1_in_int, y2_in_int, phi_jn1_rn1, alpha)
% In this function we will be computing the L1 error between the polynomial
% approximation (g(x)) and the HNA approximation (f(x)), specifically,
%
% ||g(x) - f(x)||_L1(X) = \int_{X} |g(x) - f(x)| dx, x \in X. 
%                       ~ \sum_{q = 1}^{Q} hq |g(tq) - f(tq)|. 
% 
% This may be extended to normalised error depending on the values.
%
% Inputs:
% tq - 1D integration nodes (points we are evaluating integral at)
% (xq, yq) - 2D integration nodes (points we are evaluating integral at)
% w, integration weights
% G, Screen coordinates
% L, length of screen
% k, wavenumber
% d, unit direction normal of incident wave
% n, unit normal on screen
%
% poly_coeffs, polynomial coefficients
%
% HNA_coeffs, HNA coefficients
% h_in_int, integration weights for inner integral
% (y1_in_int, y2_in_int) integration nodes for inner integral
% phi_jn1_rn1, previous iteration of phi on other screen, evaluated at the
% above integration nodes



% First computing the polynomial approximation
N_poly = length(poly_coeffs);
phi_poly = zeros(length(tq), 1);
for j = 1:length(tq)
    phi_poly(j, 1) = coeff_2_soln_midpoint_individual(poly_coeffs, L,...
        tq(j), N_poly); 
end
% computing the HNA approximation
phi_HNA = get_phi_j_r(HNA_coeffs, G, L, k, d, h_in_int, y1_in_int, ...
    y2_in_int, n, tq, xq, yq, phi_jn1_rn1, alpha);

% integration/error computation step

% err = sum(abs(phi_poly - phi_HNA).*w);

err_no_norm = sum(w.*abs(phi_poly - phi_HNA));
err_poly_only = sum(w.*abs(phi_poly));
err_HNA_only = sum(w.*abs(phi_HNA));

err_norm_poly_true = err_no_norm/err_poly_only;
err_norm_HNA_true = err_no_norm/err_HNA_only;

