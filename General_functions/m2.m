function z = m2(k, s, nq, C1, C2)
% In this function we will be computing the following kernel:
% m2(s, nq) = i H_0^(1)(k |s - nq|)/4 - sigma_1(s, t)m_1(s, t)
%
% Inputs:
% k, the wavenumber
% s, collocation point, scalar
% nq, vector of integration nodes
% C1, C2, constants needed for smoothing function
%
% Output:
% z, approximation to above integral

select = (s == nq);

dist = abs(s - nq);

z = 1i*besselh(0, k*dist)/4 ...
    - m1(k, s, nq, C1, C2).*log(k*dist); 

z(select) = 1i/4  - 2*(log(1/2) - (-psi(1)))/pi;

