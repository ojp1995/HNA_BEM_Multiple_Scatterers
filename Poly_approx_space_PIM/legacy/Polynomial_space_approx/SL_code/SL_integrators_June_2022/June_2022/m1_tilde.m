function z = m1_tilde(k,s,t,C1,C2)
% This function computes m1tilde where
% (i/4)*H_0^{(1)}(k*abs(s-t)) = log(k*abs(s-t))*m1tilde + m2
% Input parameters:
%   k is the wavenmuber
%   s is the pt of evaluation in the Hankel function (e.g. collocation pt)
%   t is the integration variable
%   C1 is the first constant for the smoothing function
%   C2 is the second constant for the smoothing function
z = (-1/(2*pi))*besselj(0,k*abs(s-t)).*smoothing_fun(k*abs(s-t),C1,C2);
end


