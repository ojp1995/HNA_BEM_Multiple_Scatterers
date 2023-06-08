function z = m1(k, s, nq, C1, C2)
% In this function we will be approximating the kernel of integral of the
% hankel function, specifically, 
%       m1(k, s, nq) = - chi(nq, C1, C2)J_0(k|s - nq|)/2pi,
% where chi(t, C1, C2) is the smoothing function and J_0 is the Bessel
% function of the first kind.
% 
% Inputs:
% k, the wavenumber
% s, the collocation points we are evluating the integral at
% nq, the quadrature nodes
% C1, C2 are constants for the smoothing function
%
% Outputs:
% z, approximation of m1

z = -smoothing_function(k*abs(s - nq), C1, C2).*besselj(0, k*abs(s - nq))/(2*pi);

end