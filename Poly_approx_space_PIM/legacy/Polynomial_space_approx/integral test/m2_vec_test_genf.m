function z = m2_vec_test_genf(k, s, t, C1, C2, f)
% This function is a test for computing the vectorised m2tilde where
% (i/4)*H_{0}^{(1)}(k*abs(s - t))*f(t) = log(k*abs(s - t))*m1tilde*f(t) + h*m2*f(t).
%
% Input paramteters:
% k is the wavenumber
% s is the point we are evaluating the integral at, the collocation point
% t is the integration variable (nodes) can be scalar or vector
% C1 and C2 are constants for the smoothing function
% f is a general function with only 1 input t. Is passed as an anon
% function  so can have more variables.

z = 1i*f(t).*besselh(0, k*abs(s - t))/4 ...
    - m1_vec_test_genf(k, s, t, C1, C2, f).*(log(k*abs(s - t)));

