function err = L1_err_poly_HNA(tq, xq, yq, w, G, L, k, d, n, poly_coeffs, HNA_coeffs,  h_in_int, y1_in_int, y2_in_int)
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



% First computing the polynomial approximation
N_poly = length(poly_coeffs);
phi_poly = zeros(length(tq), 1);
for j = 1:length(tq)
    phi_poly(j, 1) = coeff_2_soln_midpoint_individual(poly_coeffs, L,...
        tq(jh), N_poly); 
end
% computing the HNA approximation
phi_HNA = get_phi_j_r(HNA_coeffs, G, L, k, d, h_in_int, y1_in_int, ...
    y2_in_int, n, t1d, xq, yq, phi_jn1_rn1);

% integration/error computation step

err = abs(phi_poly - phi_HNA).*w;

