function I = midpoint_dphikdn_f_diff_screen(k, x1, x2, h, y1t, y2t, fnq, n)
% In this function we are integrating the normal derivate of the hankel
% function (d/dn(x)) of the usual hankel kernel with it being multiplied by
% a function. Specifically:
%
% Inputs:
% k is the wavenumber
% (x1, x2) each a vector of the collocation points
% h, midpoint weights
% (y1t, y2t), integration nodes
% fnq, function evaluated at the integration nodes
% n, normal from the screen

% Outputs:
% I , approximation to the integral.



for j = 1:length(x1)

    dist = sqrt( (x1(j) - y1t).^2 + (x2(j) - y2t).^2  );

    dist_dot_n = ( (x1(j) - y1t)*n(1) + (x2(j) - y2t)*n(2) );

    I(j, 1) = 1i*k*sum(h.*besselh(1, 1, k*dist).*dist_dot_n.*fnq./dist)/2;

%     I(j, 1) = 1i*k*h*sum(besselh(1, 1, k*dist).*dist_dot_n.*fnq./dist)/2;

    %dist_dot_n_test(j) = dist_dot_n;

    %H1_1(j) = besselh(1, 1, k*dist);

    %I_int(j) = 1i*k*h*H1_1(j).*dist_dot_n_test(j).*fnq/(2*dist);

end


