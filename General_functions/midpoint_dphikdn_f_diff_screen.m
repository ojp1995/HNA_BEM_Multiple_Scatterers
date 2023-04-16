function I = midpoint_dphikdn_f_diff_screen(k, x1, x2, h, y1t, y2t, fnq, n)
% In this function we are integrating the normal derivate of the hankel
% function (d/dn(x)) of the usual hankel kernel with it being multiplied by
% a function. Specifically:
%
% Inputs:
%
% Outputs:
%



I = 0;

for j = 1:length(x1)

    dist = sqrt( (x1(j) - y1t).^2 + (x2(j) - y2t).^2  );

    dist_dot_n = ( (x1(j) - y1t)*n(1) + (x2(j) - y2t)*n(2) );

%     dist_dot_n_test(j) = dist_dot_n;
    
    I = I + sum(1i*k*h*besselh(1, 1, k*dist).*dist_dot_n.*fnq./(2*dist));

%     H1_1(j) = besselh(1, 1, k*dist);
% 
% 
% 
%     I_int(j) = 1i*k*h*H1_1(j).*dist_dot_n_test(j).*fnq/(2*dist);

end

% I = sum(I_int);
