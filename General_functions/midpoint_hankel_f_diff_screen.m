function I = midpoint_hankel_f_diff_screen(k, x1, x2, y1t, y2t, h, fnq)
% Is this function we will be using the midpoint to approximate the
% following integral:
% 	i/4 \int_{0}^{L_2} H_{0}^{(1)}( k | \xb - \yb(t) | )f(t) d t
% \approx 
% i/4 h \sum_{q = 1}^{Q_{2}} H_{0}^{(1)}(k|\xb - \yb(n_{q})|)f(n_{q})
%
% Inputs:
    % k - wavenumber
    % (x1, x2) - this is the collocation point on $\Gamma_{1}$ we are evaluating 
    % the integral at. Both scalar
    % (y1t, y2t) - Integration nodes Vector values.
    % h is the step size/integration weights
    % fnq - function evaluated at integration nodes, is a vector value.
% Outputs:
    % I - approximation to integral 

% first computing |\xb - \yb(n_{q})|


for s = 1:length(x1)
    dist = sqrt( ( x1(s) - y1t ).^2  + (x2(s) - y2t).^2 );
    
    I(s, 1) = 1i*sum(h.*besselh(0, k*dist).*fnq)/4;

end

