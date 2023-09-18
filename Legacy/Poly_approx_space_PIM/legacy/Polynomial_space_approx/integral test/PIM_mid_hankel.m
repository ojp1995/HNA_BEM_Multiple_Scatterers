function z = PIM_mid_hankel(k, s, t_mid, t, C1, C2, f, h)
% This function will integrate the following function 
% I(f, s) = 1i/4 \int_{\Gamma} H_{0}^{1}(k|s - t|) f(t) dt  = 
% log(k*abs(s - t))*m1tilde*f(t) + h*m2tilde*f(t)
%
% For any general function f at a point s.
%
% Problem parameters
% a is the start point of the interval
% k is the wavenumber
% s is the point we are evaluating the integral at, the collocation point
% t is the integration variable (nodes) can be scalar or vector
% C1 and C2 are constants for the smoothing function
% f is a general function with only 1 input t. Is passed as an anon
% function  so can have more variables.
%
% Discretisation parameters:
% h is the step size
% N is the number of dof.

% a_disc = a + (0:N-1)*h;
% b_disc = a + (1:N)*h;




z = sum(m1_vec_test_genf(k, s, t_mid, C1, C2, f).*v1_weights_mid_vec_test(t(1:end-1),t(2:end),k,s) ...
        + h*m2_vec_test_genf(k, s, t_mid, C1, C2, f));
