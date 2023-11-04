function us = soln_in_D(X1, X2, Y1, Y2, k, w, phi)
% In this funciton we will compute the solution in the domain.
%
% Input parameter:
% (X1, X2) - coordinates in the domain
% (Y1, Y2) - quadrature nodes
% k wavenumber
% w -weights
% phi distribution on Gamma

% assumption - square domain
us = zeros(size(X1));

for q = 1:length(Y1)

    dist = sqrt( (X1 - Y1(q)).^2 + (X2 - Y2(q)).^2);

    us = us - 1i*w(q)*besselh(0, k*dist)*phi(q)/4;
end