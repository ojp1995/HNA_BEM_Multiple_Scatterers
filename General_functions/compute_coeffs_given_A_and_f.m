function v_N_Gj_r = compute_coeffs_given_A_and_f(Aj, fj_r, VHNAj)
% In this function we will compute the coefficients, translate this into an
% object that keeps the structure for the basis functions. Hopefully, we
% will also be able to export a function like object that can give the
% value for phi at any given point.
%
% Inputs:
% Aj, this is the collocation matrix, as given by Andrews code
% fj_r, this is the right hand side vector evaluated at the collocation
% points
% VHNAj, this is the structure we are using so the coefficients correspond
% to the basis functions

aj = Aj\fj_r; % computes the coefficients

v_N_Gj_r = ProjectionFunction(aj, VHNAj);




