function z = m1_vec_test(k, s, t, C1, C2)
% This function is a test for computing the vectorised m1tilde where
% (i/4)*H_{0}^{(1)}(k*abs(s - t)) = log(k*abs(s - t))*m1tilde + m2.
%
% Input paramteters:
% k is the wavenumber
% s is the point we are evaluating the integral at, the collocation point
% t is the integration variable (nodes) can be scalar or vector
% C1 and C2 are constants for the smoothing function

z = -besselj(0, k*abs(s - t)).*smoothing_fun(k*abs(s - t), C1, C2)/(2*pi);